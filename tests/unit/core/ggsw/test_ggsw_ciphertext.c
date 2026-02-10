#include "ggsw_ciphertext.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 1024
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE 1e-7

#define K_TILDEBASE 1
#define KAPPA_TILDEBASE 4
#define NLIMBS_TILDEBASE 10
#define L_TILDEBASE NLIMBS_TILDEBASE/(K_TILDEBASE+1)

//! COMMON PART (begin)

// Test ggsw_size
Test(ggsw_size, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    cr_assert(eq(i64, ggsw_size(params), NLIMBS_TILDEBASE * NLIMBSBASE));

    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

// Test ggsw_bytes
Test(ggsw_bytes, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    cr_assert(eq(i64, ggsw_bytes(params), NLIMBS_TILDEBASE * NLIMBSBASE * NBASE * 8));

    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}


//! GGSW Part (begin)

// Test ggsw_coef_number
Test(ggsw_coef_number, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    cr_assert(eq(i64, ggsw_coef_number(params), NLIMBS_TILDEBASE * NLIMBSBASE * NBASE));
    
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

// Test new_ggsw
Test(new_ggsw, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    GGSWCiphertext* ct = new_ggsw(params, NULL);

    cr_assert(eq(int, ct->mat != NULL, 1));
    cr_assert(eq(int, ct->params != NULL, 1));
    
    delete_ggsw(ct);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

// Test ggsw_Sj_Yti
Test(ggsw_Sj_Yti, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    GGSWCiphertext* ct = new_ggsw(params, NULL);

    VecBiv* vec = ggsw_Sj_Yti(params, ct->mat, 0, 2);

    // Modify the two firsts coefficients of biGLWE(-m * sk_0 / (2^kappa_tilde)^2).
    vec[0] = 1;
    vec[1] = 2;
    
    cr_assert(eq(i64, ct->mat[(K_TILDEBASE + 1)*NLIMBSBASE*NBASE],1));
    cr_assert(eq(i64, ct->mat[(K_TILDEBASE + 1)*NLIMBSBASE*NBASE + 1],2));
    
    delete_ggsw(ct);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

// Test normalize_ggsw
Test(normalize_ggsw, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    GGSWCiphertext* ct = new_ggsw(params, NULL);
    GGSWCiphertext* res = new_ggsw(params, NULL);

    normalize_ggsw(res, ct);

    delete_ggsw(ct);
    delete_ggsw(res);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}