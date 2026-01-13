#include "bivariate_polynomial.h"
#include "core/glwe/glwe.h"
#include "vec_znx_arithmetic_private.h"
#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <math.h>

#define NBASE 8
#define KBASE 8
#define KAPPABASE 32
#define NLIMBSBASE 180
#define LBASE NLIMBSBASE/(KBASE+1)

//! COMMON PART (begin)

/**
 * @brief Test poly_biv_size
 */
Test(poly_biv_size, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);

    cr_assert_eq(poly_biv_size(params), LBASE, 
                "poly_biv_size failed: got %lld, expected %lld", poly_biv_size(params), LBASE);
    
    delete_glwe_ct_params(params);
}

/**
 * @brief Test poly_biv_bytes
 */
Test(poly_univ_bytes, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);

    cr_assert_eq(poly_univ_bytes(params), NBASE*sizeof(int64_t), 
                "poly_univ_bytes failed: got %lld, expected %lld", poly_univ_bytes(params), NBASE*sizeof(int64_t));
    
    delete_glwe_ct_params(params);
}

/**
 * @brief Test poly_biv_bytes
 */
Test(poly_biv_bytes, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);

    cr_assert_eq(poly_biv_bytes(params), NBASE*LBASE*sizeof(int64_t), 
                "poly_biv_bytes failed: got %lld, expected %lld", poly_biv_bytes(params), NBASE*LBASE*sizeof(int64_t));
    
    delete_glwe_ct_params(params);
}



//! BIV POLY PART (begin) 

/**
 * @brief Tests if it returns the right size for N in [1,100]
 * @note n_limbs = (k + 1) * l
 */
Test(poly_biv_coef_number, classic_params){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE); 

    cr_assert_eq(poly_biv_coef_number(params), NBASE * LBASE, 
                "poly_biv_coef_number failed: got %lld, expected %lld", poly_biv_coef_number(params), NBASE * LBASE);

    delete_glwe_ct_params(params);
}

/**
 * @brief Test normal_bivariate_poly
 * 
 */
Test(new_normal_random_biv_poly, basic){
    MODULE* module = new_module_info(NBASE, FFT64);
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    PolyBiv* a = new_normal_random_biv_poly(module, params);

    cr_assert(eq(int, a != NULL, 1, "new_normal_random_biv returned a NULL pointer."));

    free(a);
    delete_glwe_ct_params(params);
    delete_module_info(module);
}

/**
 * @brief Test normal_bivariate_poly
 * 
 */
Test(new_normal_random_biv_poly, is_it_normal){
    MODULE* module = new_module_info(NBASE, FFT64);
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    PolyBiv* a = new_normal_random_biv_poly(module, params);

    for(int64_t p = 0 ; p < params->N ; p++){
        cr_log_info("%lld", a[(LBASE-1)*NBASE + p]);
    }

    cr_assert(1);

    free(a);
    delete_glwe_ct_params(params);
    delete_module_info(module);
}

/**
 * @brief Test add_biv_poly correctness with random normal polynomials
 */
Test(add_biv_poly, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    MODULE* module = new_module_info(params->N, FFT64);

    PolyBiv* a = new_normal_random_biv_poly(module, params);
    PolyBiv* b = new_normal_random_biv_poly(module, params);
    PolyBiv* res = malloc(poly_biv_bytes(params));

    add_biv_poly(params, res, params->N, a, params->N, b, params->N);

    for(int64_t i = 0; i < LBASE; i++){
        for(int64_t p = 0; p < NBASE; p++){
            int64_t idx = p + i * params->N;
            cr_assert(eq(dbl, res[idx], a[idx] + b[idx]),
                "add_biv_poly mismatch at index %lld: %lld + %lld = %lld, got %lld",
                (long long)idx, a[idx], b[idx], a[idx] + b[idx], res[idx]);
        }
    }
    

    free(a); free(b); free(res);
    delete_glwe_ct_params(params);
    delete_module_info(module);
}


//! BIV POLY IN DFT PART (begin)

/**
 * @brief Test add_biv_poly_dft correctness with random normal DFT polynomials
 */
Test(add_biv_poly_dft, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    MODULE* module = new_module_info(params->N, FFT64);
    
    PolyBivDFT* a = new_normal_random_biv_poly_dft(module, params);
    PolyBivDFT* b = new_normal_random_biv_poly_dft(module, params);
    PolyBivDFT* res = malloc(poly_biv_bytes(params));

    add_biv_poly_dft(params, res, params->N, a, params->N, b, params->N);

    
    for(int64_t i = 0; i < LBASE; i++){
        for(int64_t p = 0; p < NBASE; p++){
            int64_t idx = p + i * params->N;
            cr_assert(epsilon_eq(dbl, res[idx], a[idx] + b[idx], 1e-9),
                "add_biv_poly_dft mismatch at index %lld: %f + %f = %f, got %f",
                (long long)idx, a[idx], b[idx], a[idx] + b[idx], res[idx]);
        }
    }
    

    free(a); free(b); free(res);
    delete_module_info(module);
    delete_glwe_ct_params(params);
}