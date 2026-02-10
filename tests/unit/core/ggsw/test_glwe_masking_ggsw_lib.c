#include "core/ggsw/ggsw.h"
#include "rng.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#define NBASE 4
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE (KBASE + 1)*4
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -(LBASE/2 + 1)*KAPPABASE 

#define K_TILDEBASE 1
#define KAPPA_TILDEBASE 4
#define NLIMBS_TILDEBASE 10
#define L_TILDEBASE NLIMBS_TILDEBASE/(K_TILDEBASE+1)
#define SIGMA_TILDEBASE -3

//! GGSW PART (begin)

/**
 * @brief Test glwe_secret_masking_ggsw_lib. In this test, the message is drawn in Zn[X,Y], 
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y]. 
 * Moreover, the error is small enough to not affect the message. 
 */
Test(glwe_secret_masking_ggsw_lib, small_error)
{
    // The variance of the error's distribution
    double sigma = ldexp(1.0, -(LBASE/2 + 1)*KAPPABASE);
    
    MODULE* module = new_module_info_p(NBASE);
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    
    GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NULL, NBASE, KBASE);
    
    // The result bivGLWE
    VecBiv* ct_glwe = malloc(glwe_coef_number(params_glwe)*sizeof(int64_t));

    // The input message uniformly drawn in Zn[X,Y]
    PolyBiv* msg = new_uniform_random_biv_poly(module, params_glwe, LBASE/2);

    double* msg_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, msg_univ, msg);

    // The input error normaly drawn in Zn[X,Y]
    PolyBiv* err = new_normal_random_biv_poly(module, params_glwe);

    // The final phase = msg + err
    PolyBiv* phase = calloc(NBASE*LBASE, sizeof(int64_t));
    add_biv_poly(params_glwe, phase, NBASE, msg, NBASE, err, NBASE);

    // Computes the bivGLWE ciphertext
    glwe_secret_masking_ggsw_lib(module, params_glwe, ct_glwe, sk_dft, phase);

    // The computed phase in Rn[X]
    PolyBiv* phase_computed = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
    glwe_secret_demasking_ggsw_lib(params_glwe, phase_computed, sk_dft, ct_glwe);

    // The computed phase in Rn[X]
    double* phase_computed_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, phase_computed_univ, phase_computed);

    // Compare both phase in Rn[X]
    for(int64_t p = 0 ; p < NBASE ; p++){
        double diff_1 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p];
        double diff_2 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) + ceil(phase_computed_univ[p]);
        double err_length = ldexp(1.0,-(LBASE/2 + 1)*KAPPABASE);

        int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);
        cr_assert(cond, "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p, msg_univ[p], p, phase_computed_univ[p]);
    }

    free(msg); free(msg_univ); free(err); free(phase); free(ct_glwe);
    free(phase_computed); free(phase_computed_univ); 
    delete_module_info_p(module);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_secret_key_dft(sk_dft);
}

/**
 * @brief Test glwe_secret_masking_ggsw_lib. In this test, the message is drawn in Rn[X], 
 * ie. there is a 2^(-l*kappa) loss of precision for the message, going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the bivariate message. 
 */
Test(glwe_secret_masking_ggsw_lib, uniform_RnX_message)
{
    // The variance of the error's distribution
    double sigma = ldexp(1.0, -(LBASE/2 + 1)*KAPPABASE);

    MODULE* module = new_module_info_p(NBASE);
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    
    GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NULL, NBASE, KBASE);

    // The result bivGLWE
    VecBiv* ct_glwe = malloc(glwe_coef_number(params_glwe)*sizeof(int64_t));

    // The input message uniformly drawn in Rn[X]
    double* msg_univ = malloc(poly_univ_bytes(params_glwe));
    new_normal_random_vec(NBASE, msg_univ, 1, NBASE, 0.0, 0.1);

    // Computes the bivariate form
    PolyBiv* msg = malloc(poly_biv_bytes(params_glwe));
    univ_to_biv(params_glwe, msg, msg_univ);

    // The input error normaly drawn in Zn[X,Y]
    PolyBiv* err = new_normal_random_biv_poly(module, params_glwe);

    // The final phase = msg + err
    PolyBiv* phase = calloc(NBASE*LBASE, sizeof(int64_t));
    add_biv_poly(params_glwe, phase, NBASE, msg, NBASE, err, NBASE);

    // Computes the bivGLWE ciphertext
    glwe_secret_masking_ggsw_lib(module, params_glwe, ct_glwe, sk_dft, phase);

    // The computed phase in Rn[X]
    PolyBiv* phase_computed = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
    glwe_secret_demasking_ggsw_lib(params_glwe, phase_computed, sk_dft, ct_glwe);

    // The computed phase in Rn[X]
    double* phase_computed_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, phase_computed_univ, phase_computed);

    // Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p - msgComputed_p|
    // Ie, then |err_p| + 2^(-l*kappa)
    for(int64_t p = 0 ; p < NBASE ; p++){
        double diff_1 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p];
        double diff_2 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) + ceil(phase_computed_univ[p]);
        double err_length = 3*sigma + ldexp(1.0, -LBASE*KAPPABASE);

        int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);
        cr_assert(cond, "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p, msg_univ[p], p, phase_computed_univ[p]);
    }

    free(msg); free(msg_univ); free(err); free(phase); free(phase_computed); free(phase_computed_univ);
    free(ct_glwe);
    delete_module_info_p(module);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_secret_key_dft(sk_dft);
}


//! GGSW IN DFT SPACE PART (begin)

/**
 * @brief Test glwe_secret_masking_ggsw_lib. In this test, the message is drawn in Zn[X,Y], 
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y]. 
 * Moreover, the error is small enough to not affect the message. 
 */
Test(glwe_secret_masking_ggsw_lib_dft, small_error)
{
    // The variance of the error's distribution
    double sigma = ldexp(1.0, -(LBASE/2 + 1)*KAPPABASE);
    
    MODULE* module = new_module_info_p(NBASE);
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    
    GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NULL, NBASE, KBASE);
    
    // The result bivGLWE
    VecBivDFT* ct_glwe_dft = malloc(glwe_coef_number(params_glwe)*sizeof(double));

    // The input message uniformly drawn in Zn[X,Y]
    PolyBiv* msg = new_uniform_random_biv_poly(module, params_glwe, LBASE/2);

    double* msg_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, msg_univ, msg);

    // The input error normaly drawn in Zn[X,Y]
    PolyBiv* err = new_normal_random_biv_poly(module, params_glwe);

    // The final phase = msg + err
    PolyBiv* phase = calloc(NBASE*LBASE, sizeof(int64_t));
    add_biv_poly(params_glwe, phase, NBASE, msg, NBASE, err, NBASE);

    // Computes the phase in DFT space
    PolyBivDFT* phase_dft = malloc(poly_biv_bytes(params_glwe));
    vec_znx_dft_p(module, phase_dft, LBASE, phase, LBASE, NBASE);

    // Computes the bivGLWE ciphertext
    glwe_secret_masking_ggsw_lib_dft(module, params_glwe, ct_glwe_dft, sk_dft, phase_dft);

    // The computed phase in Rn[X]
    PolyBiv* phase_computed = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
    glwe_secret_demasking_ggsw_lib_dft(params_glwe, phase_computed, sk_dft, ct_glwe_dft);

    // The computed phase in Rn[X]
    double* phase_computed_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, phase_computed_univ, phase_computed);

    // Compare both phase in Rn[X]
    for(int64_t p = 0 ; p < NBASE ; p++){
        double diff_1 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p];
        double diff_2 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) + ceil(phase_computed_univ[p]);
        double err_length = 3*sigma;

        int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);
        cr_assert(cond, "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p, msg_univ[p], p, phase_computed_univ[p]);
    }

    free(msg); free(msg_univ); free(err); free(phase); free(phase_computed); free(phase_computed_univ); 
    free(ct_glwe_dft); free(phase_dft);
    delete_module_info_p(module);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_secret_key_dft(sk_dft);
}

/**
 * @brief Test glwe_secret_masking_ggsw_lib. In this test, the message is drawn in Rn[X], 
 * ie. there is a 2^(-l*kappa) loss of precision for the message, going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the bivariate message. 
 */
Test(glwe_secret_masking_ggsw_lib_dft, uniform_RnX_message)
{
    // The variance of the error's distribution
    double sigma = ldexp(1.0, -(LBASE/2 + 1)*KAPPABASE);

    MODULE* module = new_module_info_p(NBASE);
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    
    GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NULL, NBASE, KBASE);

    // The result bivGLWE
    VecBivDFT* ct_glwe_dft = malloc(glwe_coef_number(params_glwe)*sizeof(double));

    // The input message uniformly drawn in Rn[X]
    double* msg_univ = malloc(poly_univ_bytes(params_glwe));
    new_normal_random_vec(NBASE, msg_univ, 1, NBASE, 0.0, 0.1);

    // Computes the bivariate form
    PolyBiv* msg = malloc(poly_biv_bytes(params_glwe));
    univ_to_biv(params_glwe, msg, msg_univ);

    // The input error normaly drawn in Zn[X,Y]
    PolyBiv* err = new_normal_random_biv_poly(module, params_glwe);

    // The final phase = msg + err
    PolyBiv* phase = calloc(NBASE*LBASE, sizeof(int64_t));
    add_biv_poly(params_glwe, phase, NBASE, msg, NBASE, err, NBASE);

    // Computes the phase in DFT space
    PolyBivDFT* phase_dft = malloc(poly_biv_bytes(params_glwe));
    vec_znx_dft_p(module, phase_dft, LBASE, phase, LBASE, NBASE);

    // Computes the bivGLWE ciphertext
    glwe_secret_masking_ggsw_lib_dft(module, params_glwe, ct_glwe_dft, sk_dft, phase_dft);

    // The computed phase in Rn[X]
    PolyBiv* phase_computed = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
    glwe_secret_demasking_ggsw_lib_dft(params_glwe, phase_computed, sk_dft, ct_glwe_dft);

    // The computed phase in Rn[X]
    double* phase_computed_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, phase_computed_univ, phase_computed);

    // Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p - msgComputed_p|
    // Ie, then |err_p| + 2^(-l*kappa)
    for(int64_t p = 0 ; p < NBASE ; p++){
        double diff_1 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p];
        double diff_2 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) + ceil(phase_computed_univ[p]);
        double err_length = 3*sigma + ldexp(1.0,-LBASE*KAPPABASE);

        int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);
        cr_assert(cond, "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p, msg_univ[p], p, phase_computed_univ[p]);
    }
    
    free(msg); free(msg_univ); free(err); free(phase); free(phase_computed); free(phase_computed_univ);
    free(ct_glwe_dft); free(phase_dft);
    delete_module_info_p(module);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_secret_key_dft(sk_dft);
}