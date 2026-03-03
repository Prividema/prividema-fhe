#include "core/glwe/glwe.h"
#include "common/spqlios_alias.h"
#include "rng.h"

#include <stdio.h>
#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 2
#define KBASE 1
#define KAPPABASE 4
#define NLIMBSBASE 4
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -12

void printf_poly_biv(PolyBiv* pol, int64_t pol_sl){
    printf("\nPol in Zn[X,Y]: ");
    for(int64_t i = 0 ; i < LBASE ; i++)
    {
        printf("\nY^%ld : ", i);
        for(int64_t p = 0 ; p < NBASE ; p++){
            if(pol[i*pol_sl + p] < 0)
                printf("%ld X^%ld ", pol[i*pol_sl + p], p);
            else
                printf(" %ld X^%ld ", pol[i*pol_sl + p], p);
        }
    }
}

void printf_vec_poly_biv(VecBiv* pols, int64_t pols_size){
    printf("\n\nVec in Zn[X,Y]: ");
    for(int64_t j = 0 ; j < pols_size ; j++){
        printf("\n%ld-th component", j);
        printf_poly_biv(pols + j*NBASE, pols_size);
        printf("\n");
    }
}

void printf_poly_univ(PolyUniv* pol){
    for(int64_t p = 0 ; p < NBASE ; p++){
        printf(" %ld X^%ld ", pol[p], p);
    }
}

void printf_secret_key(PolyUniv** sk_values){
    printf("\n\nSecret key: ");
    for(int64_t j = 0 ; j < KBASE ; j++){
        printf("\n%ld-th component\n", j);
        printf_poly_univ(sk_values[j]);
    }
}

Test(add_mult, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
    MODULE* module = new_module_info_p(NBASE);

    // The a_i's, b drawn uniformly
    VecBiv* ct = new_uniform_random_vec((KBASE + 1) * poly_biv_coef_number(params), 3);
    printf_vec_poly_biv(ct, KBASE + 1);

    // The secret key drawn uniformly
    GLWESecretKeyDFT* sk_dft = new_uniform_glwe_secret_key_dft(NBASE, KBASE, 3);
    PolyUniv** sk_values = transform_secret_key_values_dft_to_not_dft(sk_dft->values, NBASE, KBASE);
    printf_secret_key(sk_values);

    PolyBiv* res = calloc(NBASE * LBASE , sizeof(int64_t));
    add_mult(params, module, res, ct, sk_dft);

    printf("\n\nResult pol :");
    printf_poly_biv(res, NBASE);
    printf("\n\n");

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