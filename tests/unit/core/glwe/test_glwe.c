#include "core/glwe/glwe.h"
#include "common/spqlios_alias.h"
#include "rng.h"

#include <stdio.h>
#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 8
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -12

printf_poly_biv(PolyBiv* pol){
    printf("\n %s in Zn[X,Y]: ");
    for(int64_t i = 0 ; i < LBASE ; i++)
    {
        printf("\nY^%ld  ", i);
        for(int64_t p = 0 ; p < NBASE ; p++){
            printf(" %ld X^%ld ", pol[i*NBASE + p], p);
        }
    }
}

printf_vec_poly_biv(VecBiv* pol, int64_t pol_size){
    printf("\n %s in Zn[X,Y]: ");
    for(int64_t i = 0 ; i < pol_size ; i++){
        printf("\n\n %ld-ème component\n", i);
        for(int64_t i = 0 ; i < LBASE ; i++)
        {
            printf("\nY^%ld  ", i);
            for(int64_t p = 0 ; p < NBASE ; p++){
                printf(" %ld X^%ld ", pol[i*NBASE + p], p);
            }
        }
    }
}

printf_vec_poly_univ(VecBiv* pol, int64_t pol_size){
    printf("\n %s in Zn[X,Y]: ");
    for(int64_t i = 0 ; i < pol_size ; i++){
        printf("\n\n %ld-ème component\n", i);
        for(int64_t i = 0 ; i < LBASE ; i++)
        {
            printf("\nY^%ld  ", i);
            for(int64_t p = 0 ; p < NBASE ; p++){
                printf(" %ld X^%ld ", pol[i*NBASE + p], p);
            }
        }
    }
}

Test(add_mult, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
    MODULE* module = new_module_info_p(NBASE);

    // The a_i's, b and the secret drawn uniformly
    VecBiv* ct = new_uniform_random_vec((KBASE + 1) * poly_biv_coef_number(params), 3);
    GLWESecretKeyDFT* sk_dft = new_uniform_glwe_secret_key_dft(NBASE, KBASE, 3);

    PolyBiv* res = calloc(NBASE * LBASE , sizeof(int64_t));
    add_mult(params, module, res, ct, sk_dft);

    printf_biv(res);

    free(ct); free(res);
    delete_module_info_p(module);
    delete_glwe_ct_params(params);
    delete_glwe_secret_key_dft(sk_dft);    
}

/**
 * @brief Test glwe_encrypt_priv.
 
Test(glwe_secret_masking, basic)
{
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
    MODULE* module = new_module_info_p(NBASE);

    GLWECiphertext* ct = new_glwe(params);
    GLWESecretKeyDFT* sk_dft = new_uniform_glwe_secret_key_dft(NBASE, KBASE, 3);

    // The error drawn in Zn[X,Y]
    PolyBiv* input_phase = new_normal_random_biv_poly(module, params);

    // The input phase drawn in Zn[X,Y]
    PolyBiv* phase_input = new_normal_random_biv_poly(module, params);
    glwe_secret_masking(ct, sk_dft, phase_input);

    // The input phase in Rn[X]
    double* phase_input_univ = malloc(poly_biv_bytes(params));
    biv_to_univ(params, phase_input_univ, phase_input);

    // The computed phase in Rn[X]
    double* phase_computed_univ = malloc(poly_univ_bytes(params));
    glwe_secret_demasking(phase_computed_univ, sk_dft, ct);

    // The computed phase in Zn[X,Y]
    int64_t* phase_computed = malloc(poly_biv_bytes(params));
    univ_to_biv(params, phase_computed, phase_computed_univ);
    
    // Compare both phase in Rn[X]
    for(int64_t p = 0 ; p < NBASE ; p++){
        cr_assert(epsilon_eq(dbl, phase_input_univ[p] - floor(phase_input_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]), 0.0, ldexp(1.0,-(LBASE-1)*KAPPABASE)), 
        "Equality failed at p = %ld with phase_input_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p, phase_input_univ[p], p, phase_computed_univ[p]);
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

    free(phase_input); free(phase_input_univ); free(phase_computed_univ); free(phase_computed);
    delete_glwe(ct);
    delete_module_info_p(module);
    delete_glwe_ct_params(params);
    delete_glwe_secret_key_dft(sk_dft);
}*/