#include "core/ggsw/ggsw.h"
#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#define NBASE 4
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE 0

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

Test(glwe_secret_masking_ggsw_lib, basic){
    MODULE* module = new_module_info_p(NBASE);
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
    
    GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NULL, NBASE, KBASE);
    
    VecBiv* res = malloc(glwe_coef_number(params_glwe)*sizeof(int64_t));
    PolyBiv* phase = new_normal_random_biv_poly(module, params_glwe);

    glwe_secret_masking_ggsw_lib(module, params_glwe, res, sk_dft, phase);

    free(res); free(phase);
    delete_module_info_p(module);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_secret_key_dft(sk_dft);
}

Test(glwe_secret_demasking_ggsw_lib, basic){
    MODULE* module = new_module_info_p(NBASE);
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
    
    GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NULL, NBASE, KBASE);
    
    // The result bivGLWE
    VecBiv* ct_glwe = malloc(glwe_coef_number(params_glwe)*sizeof(int64_t));

    // The input phase drawn in Zn[X,Y]
    PolyBiv* phase_input = new_normal_random_biv_poly(module, params_glwe);
    glwe_secret_masking_ggsw_lib(module, params_glwe, ct_glwe, sk_dft, phase_input);

    // The input phase in Rn[X]
    double* phase_input_univ = malloc(poly_biv_bytes(params_glwe));
    biv_to_univ(params_glwe, phase_input_univ, phase_input);

    // The computed phase in Rn[X]
    double* phase_computed_univ = malloc(poly_univ_bytes(params_glwe));
    glwe_secret_demasking_ggsw_lib(params_glwe, phase_computed_univ, sk_dft, ct_glwe);

    // The computed phase in Zn[X,Y]
    int64_t* phase_computed = malloc(poly_biv_bytes(params_glwe));
    univ_to_biv(params_glwe, phase_computed, phase_computed_univ);
    
    // Compare both phase in Rn[X]
    for(int64_t p = 0 ; p < NBASE ; p++){
        cr_assert(epsilon_eq(dbl, phase_input_univ[p] - floor(phase_input_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]), 0.0, ldexp(1.0,-(LBASE-1)*KAPPABASE)));
    }


    printf("\n phase_input in Zn[X,Y]: ");
    for(int64_t i = 0 ; i < LBASE ; i++)
    {
        printf("\nY^%ld  ", i);
        for(int64_t p = 0 ; p < NBASE ; p++){
            printf(" %ld X^%ld ", phase_input[i*NBASE + p], p);
        }
    }
    
    printf("\n phase_input in Rn[X]: ");

    for(int64_t p = 0 ; p < NBASE ; p++){
        printf(" %lf X^%ld ", phase_input_univ[p], p);
    }

    printf("\n phase_computed in Zn[X,Y]: ");
    for(int64_t i = 0 ; i < LBASE ; i++)
    {
        printf("\nY^%ld  ", i);
        for(int64_t p = 0 ; p < NBASE ; p++){
            printf(" %ld X^%ld ", phase_computed[i*NBASE + p], p);
        }
    }
    
    printf("\n phase_computed in Rn[X]: ");

    for(int64_t p = 0 ; p < NBASE ; p++){
        printf(" %lf X^%ld ", phase_computed_univ[p], p);
    }
    

    free(ct_glwe); free(phase_input); free(phase_input_univ); free(phase_computed_univ);
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