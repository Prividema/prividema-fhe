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
#define NLIMBSBASE (KBASE + 1)*2
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -(LBASE/2 + 1)*KAPPABASE 

#define K_TILDEBASE 1
#define KAPPA_TILDEBASE 4
#define NLIMBS_TILDEBASE (K_TILDEBASE + 1)*2
#define L_TILDEBASE NLIMBS_TILDEBASE/(K_TILDEBASE+1)
#define SIGMA_TILDEBASE -3


//! GGSW PART (begin)

/** The test is done without error, it is a proof of concept*/
Test(ggsw_external_product, without_error)
{    
    // GLWE and GGSW parameters. This set of GLWE parameters is for GGSW ciphertext
    double sigma = 0;
    double sigma_tilde = 0;

    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    GLWECtParams* params_glwe_tilde = new_glwe_ct_params(NBASE, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE, SIGMA_TILDEBASE);
    GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    GGSWCiphertext* ct_ggsw = new_ggsw(params_ggsw, NULL);
    GLWECiphertext* ct_glwe_tilde = new_glwe(params_glwe_tilde);
    GLWECiphertext* res = new_glwe(params_glwe);
    
    // Draws respectively uniformly the secret key in DFT space, a Zn[X] polynomial and a Zn[X,Y] polynomial
    GGSWSecretKey* sk_ggsw = new_ggsw_secret_key(NULL, NBASE, KBASE);
    sk_ggsw->values[0][0] = 1;
    sk_ggsw->values[0][1] = 0;
    GGSWSecretKeyDFT* sk_ggsw_dft = transform_ggsw_secret_key_not_dft_to_dft(module, sk_ggsw);
    GLWESecretKeyDFT* sk_glwe_dft = transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(sk_ggsw_dft);
    
    // Draws uniformly both messages
    PolyUniv* u_univ = new_uniform_random_vec(NBASE, KAPPABASE);
    PolyBiv* m = new_uniform_random_biv_poly(module, params_glwe_tilde, 1);
    
    //! Computation with function
    // Computes ct_glwe_tilde, a bivGLWE(m) using the base 2-Kappa_tilde
    glwe_secret_masking(module, ct_glwe_tilde, sk_glwe_dft, m);

    // Computes ct_ggsw, a bivGGSW(u) using the base-2Kappa
    ggsw_secret_encrypt(module, params_ggsw, ct_ggsw, sk_ggsw_dft, u_univ); 
    
    // Computes the external product of ct_glwe_tilde and ct_ggsw
    // It should result in a bivGLWE(u*m) using the base-2Kappa decomposition
    ggsw_external_product(module, res, ct_glwe_tilde, ct_ggsw);
    normalize_glwe(module, res, res);

    // Computes the result phase = u*m + err
    PolyBiv* phase_computed = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
    glwe_secret_demasking(module, phase_computed, sk_glwe_dft, res);

    // The computed phase = u*m + err in Rn[X]
    PolyUnivRnX* um_univ_computed = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, um_univ_computed, phase_computed);

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

    //! Tests equality 
    for(int64_t p = 0 ; p < NBASE ; p++)
    {
        cr_assert(eq(dbl, um_univ_computed[p], um_univ[p]), 
        "Equality failed with um_univ_computed[%ld] = %e um_univ[%ld] = %e", p, um_univ_computed[p], p, um_univ[p]);
    }

    free(m); free(u_univ); free(u_univ_dft); free(phase_computed);
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


//! GGSW PART in DFT space (begin)

/** The test is done without error, it is a proof of concept*/
Test(ggsw_external_product_dft, without_error)
{    
    // GLWE and GGSW parameters. This set of GLWE parameters is for GGSW ciphertext
    double sigma = 0;
    double sigma_tilde = 0;

    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    GLWECtParams* params_glwe_tilde = new_glwe_ct_params(NBASE, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE, SIGMA_TILDEBASE);
    GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info_p(NBASE);

    GGSWCiphertextDFT* ct_ggsw_dft = new_ggsw_dft(params_ggsw, NULL);
    GLWECiphertextDFT* ct_glwe_tilde_dft = new_glwe_dft(params_glwe_tilde);
    GLWECiphertextDFT* res_dft = new_glwe_dft(params_glwe);

    // Computes each ciphertext out of DFT space
    GLWECiphertext* res = new_glwe(params_glwe);
    
    // Draws respectively uniformly the secret key in DFT space, a Zn[X] polynomial and a Zn[X,Y] polynomial
    GGSWSecretKey* sk_ggsw = new_ggsw_secret_key(NULL, NBASE, KBASE);
    sk_ggsw->values[0][0] = 1;
    sk_ggsw->values[0][1] = 0;
    GGSWSecretKeyDFT* sk_ggsw_dft = transform_ggsw_secret_key_not_dft_to_dft(module, sk_ggsw);
    GLWESecretKeyDFT* sk_glwe_dft = transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(sk_ggsw_dft);
    
    // Draws uniformly both messages
    PolyUniv* u_univ = new_uniform_random_vec(NBASE, KAPPABASE);
    PolyBiv* m = new_uniform_random_biv_poly(module, params_glwe_tilde, 1);
    PolyBivDFT* m_dft = malloc(poly_biv_bytes(params_glwe_tilde));
    vec_znx_dft_p(module, m_dft, LBASE, m, LBASE, NBASE);

    //! Computation with function
    // Computes ct_glwe_tilde, a bivGLWE(m) using the base 2-Kappa_tilde
    glwe_secret_masking_dft(module, ct_glwe_tilde_dft, sk_glwe_dft, m_dft);

    // Computes ct_ggsw, a bivGGSW(u) using the base-2Kappa
    ggsw_secret_encrypt_dft(module, params_ggsw, ct_ggsw_dft, sk_ggsw_dft, u_univ); 
    
    // Computes the external product of ct_glwe_tilde and ct_ggsw
    // It should result in a bivGLWE(u*m) using the base-2Kappa decomposition
    ggsw_external_product_dft(module, res_dft, ct_glwe_tilde_dft, ct_ggsw_dft);

    // res out of DFT space
    vec_znx_idft_p(module, res->vec, glwe_size(params_glwe), res_dft->pvec, glwe_size(params_glwe));
    normalize_glwe(module, res, res);

    // Computes the result phase = u*m + err
    PolyBiv* phase_computed = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
    glwe_secret_demasking_dft(module, phase_computed, sk_glwe_dft, res_dft);

    // The computed phase = u*m + err in Rn[X]
    PolyUnivRnX* um_univ_computed = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, um_univ_computed, phase_computed);

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

    //! Tests equality 
    for(int64_t p = 0 ; p < NBASE ; p++)
    {
        cr_assert(eq(dbl, um_univ_computed[p], um_univ[p]), 
        "Equality failed with um_univ_computed[%ld] = %e um_univ[%ld] = %e", p, um_univ_computed[p], p, um_univ[p]);
    }

    free(m); free(m_dft); free(m_univ); free(u_univ); free(u_univ_dft); free(phase_computed);
    free(um); free(um_univ); free(um_dft); free(um_univ_computed);

    delete_glwe(res);
    delete_glwe_dft(res_dft);
    delete_glwe_dft(ct_glwe_tilde_dft);
    delete_ggsw_dft(ct_ggsw_dft);

    delete_ggsw_secret_key(sk_ggsw);
    delete_glwe_secret_key_dft(sk_glwe_dft);
    delete_ggsw_secret_key_dft(sk_ggsw_dft);
    
    delete_glwe_ct_params(params_glwe);
    delete_glwe_ct_params(params_glwe_tilde);
    delete_ggsw_ct_params(params_ggsw);
    delete_module_info_p(module);
}
