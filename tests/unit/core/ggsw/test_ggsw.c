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

// GGSW PART (begin)

Test(add_error, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
    PolyBiv* res = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
    PolyBiv* phase = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));

    add_error(params_glwe, res, phase);
    
    delete_glwe_ct_params(params_glwe);
    free(res); free(phase);
}

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
        cr_assert(epsilon_eq(dbl, msg_univ[p] - floor(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]), 0.0, ldexp(1.0,-(LBASE/2)*KAPPABASE)), 
        "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p, msg_univ[p], p, phase_computed_univ[p]);
    }

    free(msg); free(msg_univ); free(err); free(phase); free(ct_glwe);
    free(phase_computed); free(phase_computed_univ); 
    delete_module_info_p(module);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_secret_key_dft(sk_dft);
}

/**
 * @brief Test glwe_secret_masking_ggsw_lib. In this test, the message is drawn in Zn[X,Y], 
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y].
 * Moreover, the error is to big and thus affects the message. 
 */
Test(glwe_secret_masking_ggsw_lib, big_error)
{
    // The variance of the error's distribution
    double sigma = ldexp(1.0, -(LBASE/2)*KAPPABASE);

    MODULE* module = new_module_info_p(NBASE);
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    
    GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NULL, NBASE, KBASE);
    
    // The result bivGLWE
    VecBiv* ct_glwe = malloc(glwe_coef_number(params_glwe)*sizeof(int64_t));

    // The input message uniformly drawn in Zn[X,Y] with degree l/2 in Y
    PolyBiv* msg = new_uniform_random_biv_poly(module, params_glwe, LBASE/2);

    double* msg_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params_glwe, msg_univ, msg);

    // The input error = 2*2^(-lN/2)  
    PolyBiv* err = calloc(NBASE*LBASE, sizeof(int64_t));
    for(int64_t p = 0 ; p < NBASE ; p++)
        err[(LBASE/2) * NBASE + p] = 2;

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

    // Assures the error, of length (-lN/2), affects the message, of degree l/2 in Y
    for(int64_t p = 0 ; p < NBASE ; p++){
        cr_assert(epsilon_ne(dbl, msg_univ[p] - floor(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]), 0.0, ldexp(1.0,-(LBASE/2)*KAPPABASE)), 
        "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p, msg_univ[p], p, phase_computed_univ[p]);
    }

    free(msg); free(msg_univ); free(err); free(phase); free(phase_computed); free(phase_computed_univ);
    free(ct_glwe);
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
        cr_assert(epsilon_eq(dbl, msg_univ[p] - floor(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]), 0.0, 
                  ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE)), 
        "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p, msg_univ[p], p, phase_computed_univ[p]);
    }

    free(msg); free(msg_univ); free(err); free(phase); free(phase_computed); free(phase_computed_univ);
    free(ct_glwe);
    delete_module_info_p(module);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_secret_key_dft(sk_dft);
}

Test(ggsw_external_product, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GLWECtParams* params_glwe_tilde = new_glwe_ct_params(NBASE, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE, SIGMA_TILDEBASE);
    GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    GLWECiphertext* res = new_glwe(params_glwe);
    GLWECiphertext* ct_glwe = new_glwe(params_glwe_tilde);
    GGSWCiphertext* ct_ggsw = new_ggsw(params_ggsw, NULL);

    ggsw_external_product(res, ct_glwe, ct_ggsw);

    delete_glwe(res);
    delete_glwe(ct_glwe);
    delete_ggsw(ct_ggsw);
    
    delete_glwe_ct_params(params_glwe);
    delete_glwe_ct_params(params_glwe_tilde);
    delete_ggsw_ct_params(params_ggsw);
}

