#include "core/glwe/glwe.h"
#include "common/spqlios_alias.h"
#include "rng.h"

#include <stdio.h>
#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 4
#define KBASE 2
#define KAPPABASE 4
#define NLIMBSBASE (KBASE + 1)*4
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -12

/**
 * @brief Test glwe_encrypt_priv.
 */
Test(glwe_secret_masking, random_bivariate_msg)
{
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
    MODULE* module = new_module_info_p(NBASE);

    GLWECiphertext* ct = new_glwe(params);
    GLWESecretKeyDFT* sk_dft = new_uniform_glwe_secret_key_dft(NBASE, KBASE, 3);

    // The input message computed in Zn[X,Y]
    PolyBiv* msg = new_uniform_random_biv_poly(module, params); 
    printf("Input : \n");
    printf_poly_biv(msg, NBASE, NBASE, LBASE);
    printf("\nend");

    double* msg_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params, msg_univ, msg);

    // The input error drawn in Zn[X,Y]
    PolyBiv* err = new_normal_random_biv_poly(module, params);
    printf("err : \n");
    printf_poly_biv(err, NBASE, NBASE, LBASE);
    printf("\nend\n");

    // The final phase
    PolyBiv* phase = calloc(NBASE*LBASE, sizeof(int64_t));
    add_biv_poly(params, phase, NBASE, msg, NBASE, err, NBASE);

    glwe_secret_masking(ct, sk_dft, phase);

    // The computed phase in Rn[X]
    PolyBiv* phase_computed = calloc(poly_biv_coef_number(params), sizeof(int64_t));
    glwe_secret_demasking(phase_computed, sk_dft, ct);
    double* phase_computed_univ = calloc(NBASE, sizeof(double));
    biv_to_univ(params, phase_computed_univ, phase_computed);

    printf("commputed : \n");
    printf_poly_biv(phase_computed, NBASE, NBASE, LBASE);
    printf("\nend\n");

    // Compare both phase in Rn[X]
    for(int64_t p = 0 ; p < NBASE ; p++){
        cr_assert(epsilon_eq(dbl, msg_univ[p] - floor(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]), 0.0, ldexp(1.0,SIGMABASE+4)), 
        "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p, msg_univ[p], p, phase_computed_univ[p]);
    }

    free(msg_univ); free(err); free(phase_computed); free(phase_computed_univ);
    delete_glwe(ct);
    delete_module_info_p(module);
    delete_glwe_ct_params(params);
    delete_glwe_secret_key_dft(sk_dft);
}