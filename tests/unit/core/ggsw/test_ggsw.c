#include "core/ggsw/ggsw.h"
#include "core/glwe/glwe.h"
#include "core/glwe/glwe_ciphertext.h"
#include "rng.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#define NBASE 2
#define KBASE 1
#define KAPPABASE 4
#define NLIMBSBASE (KBASE + 1)*1
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -(LBASE/2 + 1)*KAPPABASE 

#define K_TILDEBASE 1
#define KAPPA_TILDEBASE 4
#define NLIMBS_TILDEBASE (K_TILDEBASE + 1)*2
#define L_TILDEBASE NLIMBS_TILDEBASE/(K_TILDEBASE+1)
#define SIGMA_TILDEBASE -3

void printf_poly_biv_b(PolyBiv* pol, int64_t pol_sl, int64_t N, int64_t l){
    for(int64_t i = 1 ; i <= l ; i++)
    {
        printf(" Y^%ld : ", i);
        for(int64_t p = 0 ; p < N ; p++){
            if(pol[(i-1)*pol_sl + p] < 0)
                printf("%ld X^%ld ", pol[(i-1)*pol_sl + p], p);
            else
                printf(" %ld X^%ld ", pol[(i-1)*pol_sl + p], p);
        }
    }
}

void printf_vec_poly_univ_b(VecBiv* pols, int64_t pols_size, int64_t N){
    for(int64_t j = 0 ; j < pols_size ; j++){
        printf(" ");
        printf_poly_univ_ZnX(pols + j*N, N);
        printf(" ");
    }
}

void printf_glwe(VecBiv* ct_glwe, GLWECtParams* params_glwe){
    // GLWE parameters
    uint64_t N = params_glwe->N;
    printf_vec_poly_univ_b(ct_glwe, glwe_size(params_glwe), N);
}

void printf_ggsw_b(MatBiv* ct_ggsw, GGSWCtParams* params_ggsw){
    printf("\n");
    // GGSW parameters
    uint64_t k_tilde = params_ggsw->k_tilde;
    uint64_t l_tilde = params_ggsw->n_limbs_tilde / (params_ggsw->k_tilde + 1);

    // GLWE parameters
    uint64_t N = params_ggsw->params_glwe->N;
    uint64_t k = params_ggsw->params_glwe->k;
    uint64_t l = params_ggsw->params_glwe->n_limbs / (params_ggsw->params_glwe->k + 1);

    for(int64_t i = 1 ; i <= l_tilde ; i++){
        for(int64_t j = 0 ; j < k_tilde + 1; j++){
            printf("\n");
            printf_glwe(ct_ggsw + (i-1)*(k_tilde + 1)*(k+1)*N*l + j*(k+1)*N*l, params_ggsw->params_glwe);
        }
    }
}

//! GGSW PART (begin)

Test(vmp, basic)
{
    MODULE* module = new_module_info(NBASE, FFT64);

    uint64_t nrows = 4;
    uint64_t ncols = 2;

    MatUniv* mat = new_uniform_random_vec(nrows*ncols*NBASE, 3);
    printf("\n Mat:\n");
    for(int64_t i = 0 ; i < nrows ; i++){
        printf("\n");
        for(int64_t j = 0 ; j < ncols ; j++){
            printf_poly_univ_ZnX(mat + i*ncols*NBASE + j*NBASE, NBASE);
        }
    }
    for(int64_t p = 0 ; p < nrows*ncols*NBASE ; p++)
        printf(" %ld ", mat[p]);
    printf("\n");
    MatUnivDFT* mat_dft = calloc(nrows*ncols*NBASE, sizeof(double));
    vec_znx_dft_p(module, mat_dft, nrows*ncols, mat, nrows*ncols, NBASE);

    VecUniv* vec = calloc(nrows*NBASE, sizeof(int64_t));
    for(int64_t i = 0 ; i < nrows ; i++){
        vec[i*NBASE] = 1;
    }
    printf("\n");
    for(int64_t i = 0 ; i < nrows ; i++){
        printf_poly_univ_ZnX(vec + i*NBASE, NBASE);
    }

    VecUnivDFT* result_dft = calloc(ncols*NBASE,sizeof(double));
    vmp_apply_dft_p(module, result_dft, ncols, vec, nrows, NBASE, mat_dft, nrows, ncols);

    VecUniv* result = calloc(ncols*NBASE, sizeof(int64_t));
    vec_znx_idft_p(module, result, ncols, result_dft, ncols);
    
    printf("\n\n");
    for(int64_t i = 0 ; i < ncols ; i++){
        printf_poly_univ_ZnX(result + i*NBASE, NBASE);
    }printf("\n\n");

    delete_module_info(module);
    free(mat); free(mat_dft); free(vec); free(result); free(result_dft);
}
/**
Test(ggsw_external_product, basic)
{    
    // GLWE and GGSW parameters. This set of GLWE parameters is for GGSW ciphertext
    double sigma = ldexp(1.0, -(8*LBASE)*KAPPABASE);
    double sigma_tilde = ldexp(1.0, -(8*LBASE)*KAPPABASE);

    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    GLWECtParams* params_glwe_tilde = new_glwe_ct_params(NBASE, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE, SIGMA_TILDEBASE);
    GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info_p(NBASE);

    // The decryption of a bivGLWE(m) should give m_dec = m + err, and |m_dec - m| <= 2^(-kappa*l) + 3*sigma with a 99% chance
    double err_length = ldexp(1.0,-(LBASE/2)*KAPPABASE) + 3*sigma;
    cr_log_info("error length = %e", err_length);
    
    GGSWCiphertext* ct_ggsw = new_ggsw(params_ggsw, NULL);
    GLWECiphertext* ct_glwe_tilde = new_glwe(params_glwe_tilde);
    GLWECiphertext* res = new_glwe(params_glwe);
    
    // Draws respectively uniformly the secret key in DFT space, a Zn[X] polynomial and a Zn[X,Y] polynomial
    GGSWSecretKey* sk_ggsw = new_ggsw_secret_key(NULL, NBASE, KBASE);
    sk_ggsw->values[0][0] = 1;
    sk_ggsw->values[0][1] = 0;
    GGSWSecretKeyDFT* sk_ggsw_dft = transform_ggsw_secret_key_not_dft_to_dft(sk_ggsw);
    GLWESecretKeyDFT* sk_glwe_dft = transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(sk_ggsw_dft);
    
    // Draws uniformly both messages
    PolyUniv* u_univ = new_uniform_random_vec(NBASE, KAPPABASE);
    // PolyBiv* m = new_uniform_random_biv_poly(module, params_glwe_tilde, 1);
    PolyBiv* m = calloc(poly_biv_coef_number(params_glwe_tilde), sizeof(int64_t));
    m[0] = 1;
    
    //! Computation with function
    // Computes ct_glwe_tilde, a bivGLWE(m) using the base 2-Kappa_tilde
    // glwe_secret_masking(ct_glwe_tilde, sk_glwe_dft, m);
    ct_glwe_tilde->vec[0] = 1;
    ct_glwe_tilde->vec[4] = 1;
    ct_glwe_tilde->vec[8] = 1;
    ct_glwe_tilde->vec[12] = 1;

    // Computes ct_ggsw, a bivGGSW(u) using the base-2Kappa
    ggsw_secret_encrypt(params_ggsw, ct_ggsw, sk_ggsw_dft, u_univ); 
    
    // Computes the external product of ct_glwe_tilde and ct_ggsw
    // It should result in a bivGLWE(u*m) using the base-2Kappa decomposition
    uint64_t nrows = ct_ggsw->params->n_limbs_tilde;
    uint64_t ncols = ct_ggsw->params->params_glwe->n_limbs;
    printf("nrows : %ld; ncols : %ld\n", nrows, ncols);
    printf("ggsw_coef_number(ct_ggsw->params) : %ld\n",ggsw_coef_number(ct_ggsw->params));

    // Computes the GGSW ciphertext in DFT space
    MatBivDFT* pmat = malloc(ggsw_coef_number(ct_ggsw->params) * sizeof(double));
    vec_znx_dft_p(module, pmat, nrows*ncols, ct_ggsw->mat, nrows*ncols, NBASE);

    // The pointer to ExternalProduct(ct_glwe, ct_ggsw)
    VecBivDFT* result = malloc(glwe_coef_number(ct_ggsw->params->params_glwe)*sizeof(double));

    // Computes ExternalProduct(ct_glwe, ct_ggsw)
    vmp_apply_dft_p(module, result, ncols, ct_glwe_tilde->vec, nrows, NBASE, pmat, nrows, ncols);
    vec_znx_idft_p(module, res->vec, ncols, result, ncols);

    free(pmat); free(result);      

    normalize_glwe(res, res);

    //printf("\n\nGLWE_TILDE(m)\n");
    printf("\n\n");
    printf_glwe(ct_glwe_tilde->vec, ct_glwe_tilde->params);
    printf_ggsw(ct_ggsw->mat, ct_ggsw->params);
    printf("\n\n");
    //printf("\n\n PRODUCT RESULT\n");
    printf_glwe(res->vec, res->params);
    printf("\n\n");


    free(m); free(u_univ);

    delete_glwe(res);
    delete_glwe(ct_glwe_tilde);
    delete_ggsw(ct_ggsw);

    delete_ggsw_secret_key(sk_ggsw);
    delete_glwe_secret_key_dft(sk_glwe_dft);
    delete_ggsw_secret_key_dft(sk_ggsw_dft);
    
    delete_glwe_ct_params(params_glwe);
    delete_glwe_ct_params(params_glwe_tilde);
    delete_ggsw_ct_params(params_ggsw);
    delete_module_info_p(module);
}*/
/**
Test(ggsw_external_product, basic)
{    
    // GLWE and GGSW parameters. This set of GLWE parameters is for GGSW ciphertext
    double sigma = ldexp(1.0, -(8*LBASE)*KAPPABASE);
    double sigma_tilde = ldexp(1.0, -(8*LBASE)*KAPPABASE);

    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    GLWECtParams* params_glwe_tilde = new_glwe_ct_params(NBASE, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE, SIGMA_TILDEBASE);
    GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info_p(NBASE);

    // The decryption of a bivGLWE(m) should give m_dec = m + err, and |m_dec - m| <= 2^(-kappa*l) + 3*sigma with a 99% chance
    double err_length = ldexp(1.0,-(LBASE/2)*KAPPABASE) + 3*sigma;
    cr_log_info("error length = %e", err_length);
    
    GGSWCiphertext* ct_ggsw = new_ggsw(params_ggsw, NULL);
    GLWECiphertext* ct_glwe_tilde = new_glwe(params_glwe_tilde);
    GLWECiphertext* res = new_glwe(params_glwe);
    
    // Draws respectively uniformly the secret key in DFT space, a Zn[X] polynomial and a Zn[X,Y] polynomial
    GGSWSecretKey* sk_ggsw = new_ggsw_secret_key(NULL, NBASE, KBASE);
    sk_ggsw->values[0][0] = 1;
    sk_ggsw->values[0][1] = 0;
    GGSWSecretKeyDFT* sk_ggsw_dft = transform_ggsw_secret_key_not_dft_to_dft(sk_ggsw);
    GLWESecretKeyDFT* sk_glwe_dft = transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(sk_ggsw_dft);
    
    // Draws uniformly both messages
    PolyUniv* u_univ = new_uniform_random_vec(NBASE, KAPPABASE);
    // PolyBiv* m = new_uniform_random_biv_poly(module, params_glwe_tilde, 1);
    PolyBiv* m = calloc(poly_biv_coef_number(params_glwe_tilde), sizeof(int64_t));
    m[0] = 1;
    
    //! Computation with function
    // Computes ct_glwe_tilde, a bivGLWE(m) using the base 2-Kappa_tilde
    // glwe_secret_masking(ct_glwe_tilde, sk_glwe_dft, m);
    ct_glwe_tilde->vec[0] = 1;
    ct_glwe_tilde->vec[2] = 1;
    ct_glwe_tilde->vec[4] = 1;
    ct_glwe_tilde->vec[6] = 1;

    // Computes ct_ggsw, a bivGGSW(u) using the base-2Kappa
    ggsw_secret_encrypt(params_ggsw, ct_ggsw, sk_ggsw_dft, u_univ); 
    
    // Computes the external product of ct_glwe_tilde and ct_ggsw
    // It should result in a bivGLWE(u*m) using the base-2Kappa decomposition
    ggsw_external_product(res, ct_glwe_tilde, ct_ggsw);
    normalize_glwe(res, res);

    // Computes the result phase = u*m + err
    PolyBiv* phase_computed = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
    glwe_secret_demasking(phase_computed, sk_glwe_dft, res);

    // The computed phase = u*m + err in Rn[X]
    PolyUnivRnX* um_univ_computed = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, um_univ_computed, phase_computed);
    // Computation with function


    //! Computation by hand
    // Computes DFT(m)
    PolyUnivDFT* u_univ_dft = malloc(NBASE*sizeof(double));
    vec_znx_dft_p(module, u_univ_dft, 1, u_univ, 1, NBASE);

    // Computes DFT(u*m)
    PolyBivDFT* um_dft = malloc(poly_biv_bytes(params_glwe));
    svp_apply_dft_p(module, um_dft, LBASE, u_univ_dft, m, LBASE, NBASE);

    // Computes u*m in ZN[X,Y]
    PolyBiv* um = malloc(poly_biv_bytes(params_glwe));
    vec_znx_idft_p(module, um, LBASE, um_dft, LBASE);
    vec_znx_normalize_base2k_p(module, KAPPA_TILDEBASE, um, LBASE, NBASE, um, LBASE, NBASE);

    // Computes u*m in Rn[X]
    double* um_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, um_univ, um);

    double* m_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, m_univ, m);

    // Computation by hand
    // printf("\n\n m in Zn[X,y] :");
    // printf_poly_biv(m, NBASE, NBASE, L_TILDEBASE);

    //printf("\n\nGLWE_TILDE(m)\n");
    printf("\n\n");
    printf_glwe(ct_glwe_tilde->vec, ct_glwe_tilde->params);
    // printf("\n\n u :\n");
    // printf_poly_univ_ZnX(u_univ, NBASE);
    //printf("\n\nGGSW(u)");
    printf_ggsw(ct_ggsw->mat, ct_ggsw->params);
    printf("\n\n");
    //printf("\n\n PRODUCT RESULT\n");
    printf_glwe(res->vec, res->params);
    printf("\n\n");

    // printf("\n\n m in Zn[X,y] :");
    // printf_poly_biv(m, NBASE, NBASE, L_TILDEBASE);

    // printf("\n\n m in Rn[X]:\n");
    // printf_poly_univ_RnX(m_univ, NBASE);

    // printf("\n\n u :\n");
    // printf_poly_univ_ZnX(u_univ, NBASE);

    // printf("\n\n u*m by hand:\n");
    // printf_poly_univ_RnX(um_univ, NBASE);

    // printf("\n\n u*m computed:\n");
    // printf_poly_univ_RnX(um_univ_computed, NBASE);

    // printf("\n\n");

    free(m); free(m_univ); free(u_univ); free(u_univ_dft); free(phase_computed);
    free(um); free(um_univ); free(um_dft); free(um_univ_computed);

    delete_glwe(res);
    delete_glwe(ct_glwe_tilde);
    delete_ggsw(ct_ggsw);

    delete_ggsw_secret_key(sk_ggsw);
    delete_glwe_secret_key_dft(sk_glwe_dft);
    delete_ggsw_secret_key_dft(sk_ggsw_dft);
    
    delete_glwe_ct_params(params_glwe);
    delete_glwe_ct_params(params_glwe_tilde);
    delete_ggsw_ct_params(params_ggsw);
    delete_module_info_p(module);
}
*/