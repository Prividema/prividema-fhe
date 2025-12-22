#include "ggsw_ciphertext.h"
#include "vec_znx_arithmetic_private.h"

/**
 * @brief Return the pointer to biGLWE(-m * sk_j * Y^i).
 * 
 * @param ct A GGSW ciphertext.
 * @param i The degree in Y of the phase = -m * sk_j * Y^i.
 * @param j The j-th component of Sk.
 * 
 * @return VecBiv*
 */
VecBiv* ggsw_Sj_Yi(GGSWCiphertext* ct, int64_t i, int64_t j){
    // GLWE parameters
    int64_t N = ct->params->params->N;
    int64_t k = ct->params->params->k;
    int64_t n_limbs = ct->params->params->n_limbs;
    int64_t l = n_limbs/(k + 1);

    // GGSW parameters
    int64_t k_tilde = ct->params->k_tilde;

    return ct->ct + i*(k_tilde + 1)*n_limbs*N + j*n_limbs*N;
}

/**
 * @brief Return the pointer to biGLWE(DFT(-m * sk_j) * Y^i) in DFT space.
 * 
 * @param ct A GGSW ciphertext in DFT space.
 * @param i The degree in Y of the phase = -m * sk_j * Y^i.
 * @param j The j-th component of Sk.
 * 
 * @return VecBivDFT*
 */
VecBivDFT* ggsw_Sj_Yi_dft(GGSWPreparedCt* ct, int64_t i, int64_t j){
    // GLWE parameters
    int64_t N = ct->params->params->N;
    int64_t k = ct->params->params->k;
    int64_t n_limbs = ct->params->params->n_limbs;
    int64_t l = n_limbs/(k + 1);

    // GGSW parameters
    int64_t k_tilde = ct->params->k_tilde;

    return ct->ct + i*(k_tilde + 1)*n_limbs*N + j*n_limbs*N;
}

/**
 * @brief Creates a bivGGSW, filled with 0.
 * 
 * @param res The result GGSW ciphertext.
 * @param params The GGSW parameters.
 */
int new_ggsw(GGSWCiphertext* res, GGSWCtParams* params, MatBiv* ct
){
    int64_t limb_size = params->params->N;
    int64_t nb_rows = params->n_limbs_tilde;
    int64_t nb_cols = params->n_limbs_tilde;

    res = malloc(sizeof(GGSWCiphertext));
    if(res == NULL)
        return -1;

    res->params = params;   

    if(ct == NULL)
        res->ct = calloc(nb_rows * nb_cols * limb_size, sizeof(int64_t));
        if(res->ct == NULL)
            return -1;
    else
        res->ct = ct;
    
    return 0;
}

/**
 * @brief Delete a bivGGSW ciphertext.
 * 
 * @param ct The bivGGSW ciphertext.
 */
void delete_ggsw(GGSWCiphertext* ct
){
    delete_ggsw_ct_params(ct->params);
    free(ct);
}

void normalize_ggsw(GGSWCiphertext* res,
                    GGSWCiphertext* ct
){
    // GLWE parameters
    int64_t N = res->params->params->N;
    int64_t n_limbs = res->params->params->n_limbs;

    // GGSW parameters
    int64_t k_tilde = res->params->k_tilde;
    int64_t n_limbs_tilde = res->params->n_limbs_tilde;

    // Matrix parameters
    int64_t nb_partial = n_limbs_tilde/(k_tilde + 1);
    int64_t nb_rows_per_partial = k_tilde + 1;

    MODULE* module = new_module(N,FFT64);

    for(int64_t i = 0 ; i < nb_partial ; i++)
    {
        for(int64_t j = 0 ; j < nb_rows_per_partial ; j++)
        {
            // The pointer to biGLWE(-m * sk_j * Y^i)
            VecBiv* res_glwe = ggsw_Sj_Yi(res->ct, i, j);
            VecBiv* ct_glwe = ggsw_Sj_Yi(ct->ct, i, j);
            
            // Normalize ct
            uint8_t* tmp = vec_znx_normalize_base2k_tmp_bytes(module);
            vec_znx_normalize_base2k(module, ct->params->params->kappa, res_glwe, n_limbs, N, ct_glwe, n_limbs, N, tmp);
            
            free(tmp);
        }
    }

    free(module);
}

/**
 * @brief Adds two GGSW ciphertexts with same params and put result in res
 * 
 * @param res The result GGSW ciphertext. 
 * @param ct1 The left-hand side GGSW ciphertext.
 * @param ct2 The right-hand side GGSW ciphertext.
 */
void ggsw_add(GGSWCiphertext* res,  // result
             GGSWCiphertext* ct1,  // first operand
             GGSWCiphertext* ct2   // second operand
){
    int64_t nb_rows = res->params->n_limbs_tilde;
    int64_t nb_cols = res->params->params->n_limbs;
    int64_t N = res->params->params->N;

    for (int64_t i = 0 ; i < nb_rows ; i++){
        for (int64_t j = 0 ; j < nb_cols ; j++){
            for (int64_t k = 0 ; k < N ; k++){
                res->ct[i*N*nb_cols + j*N + k] = ct1->ct[i*N*nb_cols + j*N + k] + ct2->ct[i*N*nb_cols + j*N + k];
            } 
        }
    }
}

/**
 * @brief  Multiply a GGSW ciphertext by a constant in ZnX
 * 
 * @param res The result GGSW ciphertext.
 * @param ct1 The GGSW ciphertext.
 * @param u The polynomial in ZnX, with coefficient in [-2^(kappa-1), 2^(kappa-1)]
 */
void ggsw_const_mult(GGSWCiphertext* res,  
                     GGSWCiphertext* ct, 
                     PolyUniv* u
){
    int64_t nb_rows = res->params->n_limbs_tilde;
    int64_t nb_cols = res->params->params->n_limbs;
    int64_t N = res->params->params->N;

    MODULE* module = new_module_info(N, FFT64);
    
    // The polynomial in DFT space
    SVP_PPOL* u_prep = new_svp_ppol(module);
    svp_prepare(module, u_prep, u);

    // The ciphertext in DFT space
    MatBivDFT_* ct_dft = new_vec_znx_dft(module, nb_rows * nb_cols);
    vec_znx_dft(module, ct_dft, nb_rows * nb_cols, ct->ct, nb_rows * nb_cols, N);

    svp_apply_dft(module, ct_dft, nb_rows * nb_cols, u_prep, ct, nb_rows * nb_cols, N);

    // Go back to ZnXY
    uint8_t* temp = NULL;
    vec_znx_idft(module, (VEC_ZNX_BIG*)res->ct, nb_rows * nb_cols, ct_dft, nb_rows * nb_cols, temp);
    
    // Normalization
    for(int64_t i = 0 ; i < nb_partial ; i++)
    {
        for()
        {
            // The pointer to biGLWE(-m * sk_j * Y^i)
            vec_znx_normalize_base2k(module, ct->params->params->kappa, 
                                    );
        }
    }
}


/**
 * @brief Creates a bivGGSW in DFT space, filled with 0.
 * 
 * @param res The result bivGGSW ciphertext in DFT space.
 * @param params The GGSW parameters.
 */
int new_ggsw_prepared(GGSWPreparedCt* res, GGSWCtParams* params, MatBivDFT* ct
){
    int64_t limb_size = params->params->N;
    int64_t nb_rows = params->n_limbs_tilde;
    int64_t nb_cols = params->n_limbs_tilde;

    res = malloc(sizeof(GGSWCiphertext));
    if(res == NULL)
        return -1;

    res->params = params;   

    if(ct == NULL)
        res->ct = calloc(nb_rows * nb_cols * limb_size, sizeof(double));
        if(res->ct == NULL)
            return -1;
    else
        res->ct = ct;
    
    return 0;
}

/**
 * @brief Delete a bivGGSW in DFT space.
 * 
 * @param ct The bivGGSW in DFT space.
 */
void delete_ggsw_prepared(GGSWPreparedCt* ct
){
    delete_ggsw_ct_params(ct->params);
    free(ct);
}

/**
 * @brief  Multiply a GGSW ciphertext by a constant in ZnX
 * 
 * @param res The result GGSW ciphertext.
 * @param ct1 The GGSW ciphertext.
 * @param u The polynomial in ZnX, with coefficient in [-2^(kappa-1), 2^(kappa-1)]
 */
void ggsw_const_mult_dft(GGSWPreparedCt* res,  
                         GGSWPreparedCt* ct, 
                         SVP_PPOL* u
){
    int64_t nb_rows = res->params->n_limbs_tilde;
    int64_t nb_cols = res->params->params->n_limbs;
    int64_t N = res->params->params->N;

    for (int64_t i = 0 ; i < nb_rows ; i++){
        for (int64_t j = 0 ; j < nb_cols ; j++){
            // Passing in DFT space
            svp_apply_dft();
            for (int64_t k = 0 ; k < N ; k++){
                res->ct[i*N*nb_cols + j*N + k] = cst * ct->ct[i*N*nb_cols + j*N + k];
            } 
        }
    }
}
