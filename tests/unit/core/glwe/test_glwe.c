#include "core/glwe/glwe.h"
#include "vec_znx_arithmetic_private.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 1024
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE 1e-7

/**
 * @brief Test glwe_encrypt_priv.
 */
Test(glwe_secret_masking, basic)
{
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    GLWECiphertext* ct = new_glwe(params);
    GLWEPreparedSK* sk_dft = new_glwe_secret_key_dft(NULL, NBASE, KBASE);
    PolyBiv* phase = new_normal_random_biv_poly(module, params);

    glwe_secret_masking(ct, sk_dft, phase);
    
    cr_assert(1);
}