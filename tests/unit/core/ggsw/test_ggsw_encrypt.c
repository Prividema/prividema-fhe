#include "core/ggsw/ggsw.h"
#include "rng.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#define NBASE 4
#define KBASE 1
#define KAPPABASE 4
#define NLIMBSBASE (KBASE + 1)*4
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -(LBASE/2 + 1)*KAPPABASE 

#define K_TILDEBASE 1
#define KAPPA_TILDEBASE 4
#define NLIMBS_TILDEBASE 10
#define L_TILDEBASE NLIMBS_TILDEBASE/(K_TILDEBASE+1)
#define SIGMA_TILDEBASE -3

Test(ggsw_secret_encrypt, does_run)
{
    // GLWE and GGSW parameters. This set of GLWE parameters is for GGSW ciphertext
    double sigma = ldexp(1.0, SIGMABASE);
    GLWECtParams* params_glwe_for_ggsw = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe_for_ggsw, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    GGSWCiphertext* ct_ggsw = new_ggsw(params_ggsw, NULL);
    GGSWSecretKeyDFT* sk_dft = new_uniform_ggsw_secret_key_dft(NBASE, KBASE, 3);
    PolyUniv* msg = new_uniform_random_vec(NBASE, 3);

    ggsw_secret_encrypt(params_ggsw, ct_ggsw, sk_dft, msg);


    free(msg);
    delete_ggsw(ct_ggsw);
    delete_ggsw_secret_key_dft(sk_dft);
    delete_glwe_ct_params(params_glwe_for_ggsw);
    delete_ggsw_ct_params(params_ggsw);
}