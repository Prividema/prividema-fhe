#include "ggsw_ciphertext.h"
#include "vec_znx_arithmetic_private.h"
#include "spqlios_alias.h"

//! GGSW Part (begin)

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
        res->mat = calloc(nb_rows * nb_cols * limb_size, sizeof(int64_t));
        if(res->mat == NULL)
            return -1;
    else
        res->mat = ct;
    
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

    return ct->mat + i*(k_tilde + 1)*n_limbs*N + j*n_limbs*N;
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

    MODULE* module = new_module_info(N,FFT64);

    for(int64_t i = 0 ; i < nb_partial ; i++)
    {
        for(int64_t j = 0 ; j < nb_rows_per_partial ; j++)
        {
            // The pointer to biGLWE(-m * sk_j * Y^i)
            VecBiv* res_glwe = ggsw_Sj_Yi(res, i, j);
            VecBiv* ct_glwe = ggsw_Sj_Yi(ct, i, j);
            
            // Normalize ct
            uint8_t* tmp = malloc(vec_znx_normalize_base2k_tmp_bytes(module));
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
                res->mat[i*N*nb_cols + j*N + k] = ct1->mat[i*N*nb_cols + j*N + k] + ct2->mat[i*N*nb_cols + j*N + k];
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
    // GGSW & GLWE params
    GGSWCtParams* params_ggsw = res->params;
    GLWECtParams* params_glwe = params_ggsw->params;

    int64_t N = res->params->params->N;
    int64_t mat_size = ggsw_size(params_ggsw);
    MODULE* module = new_module_info(N, FFT64);
    
    // The polynomial in DFT space
    PolyUnivDFT* u_prep = new_svp_ppol_p(module);
    svp_prepare_p(module, u_prep, u);

    // The ciphertext in DFT space
    MatBivDFT* ct_dft = new_vec_znx_dft_p(module, mat_size);
    vec_znx_dft_p(module, ct_dft, mat_size, ct->mat, mat_size, N);

    svp_apply_dft_p(module, ct_dft, mat_size, u_prep, ct->mat, mat_size, N);

    // Go back to ZnXY
    uint8_t* tmp = NULL;
    vec_znx_idft_p(module, res->mat, mat_size, ct_dft, mat_size, tmp);
    
    // Normalization
    for(int64_t i = 0 ; i < nb_partials(params_ggsw) ; i++)
    {
        for(int64_t j = 0 ; j < nb_rows_per_partial(params_ggsw) ; j++)
        {
            // The pointer to biGLWE(-m * sk_j * Y^i)
            VecBiv* ct_biv = ggsw_Sj_Yi(res, i, j);

            // TODO Does it works to do it inplace ?
            uint8_t* tmp = malloc(vec_znx_big_normalize_base2k_tmp_bytes(module));
            vec_znx_normalize_base2k(module, ct->params->params->kappa, 
                                     ct_biv, glwe_size(params_glwe), N,
                                     ct_biv, glwe_size(params_glwe), N,
                                     tmp
                                    );
            free(tmp);
        }
    }
}


//! GGSW DFT PART (begin)

/**
 * @brief Return the pointer to biGLWE(DFT(-m * sk_j) * Y^i) in DFT space.
 * 
 * @param ct A GGSW ciphertext in DFT space.
 * @param i The degree in Y of the phase = -m * sk_j * Y^i.
 * @param j The j-th component of Sk.
 * 
 * @return VecBivDFT*
 */
VecBivDFT* ggsw_Sj_Yi_dft(GGSWCiphertextDFT* ct_dft, int64_t i, int64_t j){
    // GLWE parameters
    int64_t N = ct_dft->params->params->N;
    int64_t k = ct_dft->params->params->k;
    int64_t n_limbs = ct_dft->params->params->n_limbs;
    int64_t l = n_limbs/(k + 1);

    // GGSW parameters
    int64_t k_tilde = ct_dft->params->k_tilde;

    return ct_dft->pmat + i*(k_tilde + 1)*n_limbs*N + j*n_limbs*N;
}




/**
 * @brief Creates a bivGGSW in DFT space, filled with 0.
 * 
 * @param res The result bivGGSW ciphertext in DFT space.
 * @param params The GGSW parameters.
 */
int new_ggsw_prepared(GGSWCiphertextDFT* res_dft, GGSWCtParams* params, MatBivDFT* ct
){
    int64_t limb_size = params->params->N;
    int64_t nb_rows = params->n_limbs_tilde;
    int64_t nb_cols = params->n_limbs_tilde;

    res_dft = malloc(sizeof(GGSWCiphertext));
    if(res_dft == NULL)
        return -1;

    res_dft->params = params;   

    if(ct == NULL)
        res_dft->pmat = calloc(nb_rows * nb_cols * limb_size, sizeof(double));
        if(res_dft->pmat == NULL)
            return -1;
    else
        res_dft->pmat = ct;
    
    return 0;
}

/**
 * @brief Delete a bivGGSW in DFT space.
 * 
 * @param ct The bivGGSW in DFT space.
 */
void delete_ggsw_prepared(GGSWCiphertextDFT* res_dft
){
    delete_ggsw_ct_params(res_dft->params);
    free(res_dft);
}

/**
 * @brief  Multiply a GGSW ciphertext by a constant in ZnX
 * 
 * @param res The result GGSW ciphertext.
 * @param ct1 The GGSW ciphertext.
 * @param u The polynomial in ZnX, with coefficient in [-2^(kappa-1), 2^(kappa-1)]
 */
/**
 * @brief  Multiply a GGSW ciphertext by a constant in ZnX
 * 
 * @param res The result GGSW ciphertext.
 * @param ct1 The GGSW ciphertext.
 * @param u The polynomial in ZnX, with coefficient in [-2^(kappa-1), 2^(kappa-1)]
 */
void ggsw_const_mult_dft(GGSWCiphertextDFT* res_dft,  
                         GGSWCiphertextDFT* ct_dft, 
                         PolyUniv* u
){
    // GGSW & GLWE params
    GGSWCtParams* params_ggsw = res_dft->params;
    GLWECtParams* params_glwe = params_ggsw->params;

    int64_t N = res_dft->params->params->N;
    int64_t mat_size = ggsw_size(params_ggsw);
    MODULE* module = new_module_info(N, FFT64);
    
    // Does u_dft = DFT(u)
    PolyUnivDFT* u_dft = new_svp_ppol_p(module);
    svp_prepare_p(module, u_dft, u);

    // Temporary GGSW ciphertext
    MatBiv* tmp_mat = NULL;
    GGSWCiphertext* tmp_ggsw_1;
    new_ggsw(tmp_ggsw_1, params_ggsw, tmp_mat);

    // Does tmp_ggsw = iDFT(ct_in_dft). Then ct_dft = DFT(u) * DFT(iDFT(ct_dft))) = DFT(u) * ct_in_dft
    uint8_t* tmp_space = NULL;
    vec_znx_idft_p(module, tmp_ggsw_1->mat, mat_size, ct_dft->pmat, mat_size, tmp_space);
    svp_apply_dft_p(module, ct_dft->pmat, mat_size, u_dft, tmp_ggsw_1->mat, mat_size, N);

    delete_ggsw(tmp_ggsw_1);

    // Does tmp_ggsw = iDFT(ct_dft) = u * iDFT(ct_dft_in)
    MatBiv* tmp = NULL;
    GGSWCiphertext* tmp_ggsw_2;
    new_ggsw(tmp_ggsw_2, params_ggsw, tmp);

    vec_znx_idft_p(module, tmp_ggsw_2->mat, mat_size, ct_dft->pmat, mat_size, tmp);
    
    // Normalization of tmp_ggsw = u * iDFT(ct_dft_in) 
    for(int64_t i = 0 ; i < nb_partials(params_ggsw) ; i++)
    {
        for(int64_t j = 0 ; j < nb_rows_per_partial(params_ggsw) ; j++)
        {
            // The pointer to biGLWE(-m * sk_j * Y^i)
            VecBivDFT* ct_biv = ggsw_Sj_Yi(tmp_ggsw_2, i, j);

            // TODO Does it works to do it inplace ?
            uint8_t* tmp = vec_znx_big_normalize_base2k_tmp_bytes(module);
            vec_znx_normalize_base2k(module, ct_dft->params->params->kappa, 
                                     ct_biv, glwe_size(params_glwe), N,
                                     ct_biv, glwe_size(params_glwe), N,
                                     tmp
                                    );
            free(tmp);
        }
    }

    // Go back to DFT space
    vec_znx_dft_p(module, res_dft->pmat, mat_size, tmp_ggsw_2, mat_size, N);

    delete_ggsw(tmp_ggsw_2);
}


//! COMMON PART (begin)

/**
 * @brief Return the size of a bivGGSW ciphertext, in DFT space & out of DFT space.
 * 
 * @param params The GGSW parameters.
 * @return int64_t 
 * 
 * @note The size of a bivGGSW ciphertext is the same in and out of DFT space.
 */
int64_t ggsw_size(GGSWCtParams* params){
    return params->n_limbs_tilde * params->params->n_limbs;
}

/**
 * @brief The number of bytes needed to store a bivGGSW ciphertext.
 * 
 * @param params The GGSW parameters.
 * @return int64_t 
 * 
 * @note The number of bytes needed to store a bivGGSW ciphertext, is the same in and out of DFT space. 
 */
int64_t ggsw_bytes(GGSWCtParams* params){
    int64_t N = params->params->N;
    return ggsw_size(params) * N * sizeof(int64_t); 
}



