#include "core/glwe/glwe.h"
#include "vec_znx_arithmetic_private.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 8
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -12

/**
 * @brief Test glwe_encrypt_priv.
 */
Test(glwe_secret_masking, basic)
{
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
    MODULE* module = new_module_info(NBASE, FFT64);

    GLWECiphertext* ct = new_glwe(params);
    GLWEPreparedSK* sk_dft = new_uniform_glwe_secret_key_dft(NBASE, KBASE, 3);

    // Draws err in Zn[X,Y]
    PolyBiv* input_phase = new_normal_random_biv_poly(module, params);

    // Computes bivGLWE(msg + err)
    glwe_secret_masking(ct, sk_dft, input_phase);
    /*
    // Computes err
    PolyBiv* computed_phase = malloc(poly_biv_bytes(params));
    glwe_secret_demasking(params, computed_phase, sk_dft, ct);

    for(int64_t i = 0 ; i < LBASE ; i++){
        for(int64_t p = 0 ; p < NBASE ; p++){
            cr_assert(eq(i64, computed_phase[i*NBASE + p], input_phase[i*NBASE + p]));
        }
    }*/
}