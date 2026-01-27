#include "ggsw_ciphertext.h"
#include "vec_znx_arithmetic_private.h"
#include "spqlios_alias.h"
#include "glwe_ciphertext.h"

//! GGSW Part (begin)

uint64_t ggsw_coef_number(GGSWCtParams* params){
uint64_t ggsw_coef_number(GGSWCtParams* params){
    return params->n_limbs_tilde * glwe_coef_number(params->params);
}

GGSWCiphertext* new_ggsw(GGSWCtParams* params, MatBiv* mat)
{
    uint64_t limb_size = params->params->N;
    uint64_t nb_rows = params->n_limbs_tilde;
    uint64_t nb_cols = params->n_limbs_tilde;
    uint64_t limb_size = params->params->N;
    uint64_t nb_rows = params->n_limbs_tilde;
    uint64_t nb_cols = params->n_limbs_tilde;

    GGSWCiphertext* ggsw_ct = malloc(sizeof(GGSWCiphertext));
    if(ggsw_ct == NULL){
        perror("Malloc failed.");
        return NULL;
    }
    
    ggsw_ct->params = params;   

    if(mat == NULL)
    {
        ggsw_ct->mat = calloc(ggsw_coef_number(params), sizeof(int64_t));
        if(ggsw_ct->mat == NULL){
            perror("Malloc failed.");
            return NULL;
        }
    }
    else
        ggsw_ct->mat = mat;
    
    return ggsw_ct;
}

void delete_ggsw(GGSWCiphertext* ct)
{
    free(ct->mat);
    free(ct);
}

VecBiv* ggsw_Sj_Yti(GGSWCiphertext* ct, int64_t j, int64_t i)
{
    // GLWE parameters
    uint64_t N = ct->params->params->N;
    uint64_t k = ct->params->params->k;
    uint64_t n_limbs = ct->params->params->n_limbs;
    uint64_t l = n_limbs/(k + 1);

    // GGSW parameters
    uint64_t k_tilde  = ct->params->k_tilde;

    return ct->mat + i*(k_tilde + 1)*n_limbs*N + j*n_limbs*N;
}

void normalize_ggsw(GGSWCiphertext* res, GGSWCiphertext* ct)
{
    // GLWE parameters
    uint64_t N = res->params->params->N;
    uint64_t n_limbs = res->params->params->n_limbs;

    // GGSW parameters
    uint64_t k_tilde  = res->params->k_tilde;
    uint64_t n_limbs_tilde = res->params->n_limbs_tilde;

    // Matrix parameters
    uint64_t nb_partial = n_limbs_tilde/(k_tilde + 1);
    uint64_t nb_rows_per_partial = k_tilde + 1;

    MODULE* module = new_module_info(N,FFT64);

    for(int64_t i = 0 ; i < nb_partial ; i++)
    {
        for(int64_t j = 0 ; j < nb_rows_per_partial ; j++)
        {
            // The pointer to biGLWE(-m * sk_j * Y^i)
            VecBiv* res_glwe = ggsw_Sj_Yti(res, j, i);
            VecBiv* ct_glwe = ggsw_Sj_Yti(ct, j, i);
            
            // Normalize ct
            vec_znx_normalize_base2k_p(module, ct->params->params->kappa, res_glwe, n_limbs, N, ct_glwe, n_limbs, N);
        }
    }
    free(module);
}

void add_ggsw(GGSWCiphertext* res,  // result
             GGSWCiphertext* ct1,   // first operand
             GGSWCiphertext* ct2)   // second operand
{
    uint64_t nb_rows = res->params->n_limbs_tilde;
    uint64_t nb_cols = res->params->params->n_limbs;
    uint64_t N = res->params->params->N;

    for (int64_t i = 0 ; i < nb_rows ; i++)
        for (int64_t j = 0 ; j < nb_cols ; j++)
            for (uint64_t k = 0 ; k < N ; k++)
                res->mat[i*N*nb_cols + j*N + k] = ct1->mat[i*N*nb_cols + j*N + k] + ct2->mat[i*N*nb_cols + j*N + k];
}

void const_mult_ggsw(GGSWCiphertext* res,  
                     GGSWCiphertext* ct, 
                     PolyUniv* u)
{
    // GGSW & GLWE params
    GGSWCtParams* params_ggsw = res->params;
    GLWECtParams* params_glwe = params_ggsw->params;

    uint64_t N = res->params->params->N;
    int64_t mat_size = ggsw_size(params_ggsw);
    MODULE* module = new_module_info(N, FFT64);
    
    // The polynomial in DFT space
    PolyUnivDFT* u_prep = new_svp_ppol_p(module);
    svp_prepare_p(module, u_prep, u);

    // The ciphertext in DFT space
    MatBivDFT* ct_dft = new_vec_znx_dft_p(module, mat_size);
    vec_znx_dft_p(module, ct_dft, mat_size, ct->mat, mat_size, N);

    svp_apply_dft_p(module, ct_dft, mat_size, u_prep, ct->mat, mat_size, N);

    // Go back to Zn[X,Y]
    vec_znx_idft_p(module, res->mat, mat_size, ct_dft, mat_size);
    
    // Normalization
    for(int64_t i = 0 ; i < nb_partials(params_ggsw) ; i++)
    {
        for(int64_t j = 0 ; j < nb_rows_per_partial(params_ggsw) ; j++)
        {
            // The pointer to biGLWE(-m * sk_j * Y^i)
            VecBiv* ct_biv = ggsw_Sj_Yti(res, j, i);

            // TODO Does it works to do it inplace ?
            vec_znx_normalize_base2k_p(module, ct->params->params->kappa, 
                                     ct_biv, glwe_size(params_glwe), N,
                                     ct_biv, glwe_size(params_glwe), N
            );
        }
    }
}


//! GGSW DFT PART (begin)

uint64_t ggsw_coef_number_dft(GGSWCtParams* params){
uint64_t ggsw_coef_number_dft(GGSWCtParams* params){
    return (params->n_limbs_tilde * glwe_coef_number(params->params))/2;
}

GGSWCiphertextDFT* new_ggsw_prepared(GGSWCtParams* params, MatBivDFT* pmat)
{
    GGSWCiphertextDFT* ggsw_ct_dft = malloc(sizeof(GGSWCiphertext));
    if(ggsw_ct_dft == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    ggsw_ct_dft->params = params;   

    if(pmat == NULL)
    {
        ggsw_ct_dft->pmat = calloc(2 * ggsw_coef_number_dft(params), sizeof(double));
        if(ggsw_ct_dft->pmat == NULL){
            perror("Malloc failed.");
            return NULL;
        }
    }
    else
        ggsw_ct_dft->pmat = pmat;
    
    return 0;
}

void delete_ggsw_dft(GGSWCiphertextDFT* res_dft)
{
    free(res_dft->pmat);
    free(res_dft);
}

VecBivDFT* ggsw_Sj_Yti_dft(GGSWCiphertextDFT* ct_dft, int64_t j, int64_t i)
{
    // GLWE parameters
    uint64_t N = ct_dft->params->params->N;
    uint64_t k = ct_dft->params->params->k;
    uint64_t n_limbs = ct_dft->params->params->n_limbs;
    uint64_t l = n_limbs/(k + 1);

    // GGSW parameters
    uint64_t k_tilde  = ct_dft->params->k_tilde;

    return ct_dft->pmat + i*(k_tilde + 1)*n_limbs*N + j*n_limbs*N;
}

void const_mult_ggsw_dft(GGSWCiphertextDFT* res_dft,  
                         GGSWCiphertextDFT* ct_dft, 
                         PolyUniv* u)
{
    // GGSW & GLWE params
    GGSWCtParams* params_ggsw = res_dft->params;
    GLWECtParams* params_glwe = params_ggsw->params;

    uint64_t N = res_dft->params->params->N;
    int64_t mat_size = ggsw_size(params_ggsw);
    MODULE* module = new_module_info(N, FFT64);
    
    // Does u_dft = DFT(u)
    PolyUnivDFT* u_dft = new_svp_ppol_p(module);
    svp_prepare_p(module, u_dft, u);

    // Temporary GGSW ciphertext
    GGSWCiphertext* tmp_ggsw_1 = new_ggsw(params_ggsw, NULL);

    // Does tmp_ggsw = iDFT(ct_in_dft). Then ct_dft = DFT(u) * DFT(iDFT(ct_dft))) = DFT(u) * ct_in_dft
    vec_znx_idft_p(module, tmp_ggsw_1->mat, mat_size, ct_dft->pmat, mat_size);
    svp_apply_dft_p(module, ct_dft->pmat, mat_size, u_dft, tmp_ggsw_1->mat, mat_size, N);

    delete_ggsw(tmp_ggsw_1);

    // Does tmp_ggsw = iDFT(ct_dft) = u * iDFT(ct_dft_in)
    GGSWCiphertext* tmp_ggsw_2 = new_ggsw(params_ggsw, NULL);

    vec_znx_idft_p(module, tmp_ggsw_2->mat, mat_size, ct_dft->pmat, mat_size);
    
    // Normalization of tmp_ggsw = u * iDFT(ct_dft_in) 
    for(int64_t i = 0 ; i < nb_partials(params_ggsw) ; i++)
    {
        for(int64_t j = 0 ; j < nb_rows_per_partial(params_ggsw) ; j++)
        {
            // The pointer to biGLWE(-m * sk_j * Y^i)
            VecBiv* ct_biv = ggsw_Sj_Yti(tmp_ggsw_2, j, i);

            // TODO Does it works to do it inplace ?
            vec_znx_normalize_base2k_p(module, ct_dft->params->params->kappa, 
                                     ct_biv, glwe_size(params_glwe), N,
                                     ct_biv, glwe_size(params_glwe), N
            );
        }
    }

    // Go back to DFT space
    vec_znx_dft_p(module, res_dft->pmat, mat_size, tmp_ggsw_2->mat, mat_size, N);

    delete_ggsw(tmp_ggsw_2);
}


//! COMMON PART (begin)

uint64_t ggsw_size(GGSWCtParams* params){
uint64_t ggsw_size(GGSWCtParams* params){
    return params->n_limbs_tilde * params->params->n_limbs;
}

uint64_t ggsw_bytes(GGSWCtParams* params){
    int64_t N = params->params->N;
    return ggsw_size(params) * N * sizeof(int64_t); 
}
