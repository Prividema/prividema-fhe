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

Test(add_error, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
    PolyBiv* res = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
    PolyBiv* phase = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));

    add_error(params_glwe, res, phase);
    
    delete_glwe_ct_params(params_glwe);
    free(res); free(phase);
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

