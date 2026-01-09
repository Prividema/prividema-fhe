#include "bivariate_polynomial.h"
#include "core/glwe/glwe.h"
#include "vec_znx_arithmetic_private.h"
#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define NBASE 1024
#define KBASE 8
#define KAPPABASE 32
#define NLIMBSBASE 160

//! BIV POLY PART (begin) 

/**
 * @brief Tests if it returns the right size for N in [1,100]
 * @note n_limbs = (k + 1) * l
 */
Test(poly_biv_coef_number, classic_params){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE); 

    cr_assert_eq(poly_biv_coef_number(params), NLIMBSBASE/(KBASE + 1), "poly_biv_size failed: got %lld, expected 4", poly_biv_size(params));

    delete_glwe_ct_params(params);
}


int normal_random_biv_poly(MODULE* module, 
                           GLWECtParams*  params, 
                           PolyBiv* res
);

/**
 * @brief Test normal_bivariate_poly
 * 
 */
Test(normal_bivariate_poly, basic){
    MODULE* module = new_module_info(NBASE, FFT64);
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    PolyBiv* pol = normal_bivariate_poly(module, params, );
    cr_assert_eq(poly_biv_size(params), n_limbs/2, "poly_biv_size failed: got %lld, expected 4", poly_biv_size(params));
}

/**
 * @brief Test add_biv_poly correctness with random normal polynomials
 */
Test(add_biv_poly, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    int64_t size = poly_biv_size(params);
    int64_t total = size * params->N;

    PolyBiv* a = malloc(sizeof(PolyBiv) * total);
    PolyBiv* b = malloc(sizeof(PolyBiv) * total);
    PolyBiv* res = malloc(sizeof(PolyBiv) * total);

    MODULE* module = new_module_info(params->N, FFT64);

    normal_random_biv_poly(module, params, a);
    normal_random_biv_poly(module, params, b);

    add_biv_poly(params, res, params->N, a, params->N, b, params->N);

    for(int64_t i = 0; i < size; i++){
        for(int64_t j = 0; j < params->N; j++){
            int64_t idx = j + i * params->N;
            cr_assert_eq(res[idx], a[idx] + b[idx],
                "add_biv_poly mismatch at index %lld: %lld + %lld = %lld, got %lld",
                (long long)idx, (long long)a[idx], (long long)b[idx], (long long)(a[idx] + b[idx]), (long long)res[idx]);
        }
    }

    free(a); free(b); free(res);
    delete_glwe_ct_params(params);
    delete_module_info(module);
}



//! BIV POLY IN DFT PART (begin)

int64_t poly_biv_coef_number_dft(GLWECtParams* params);

void add_biv_poly_dft(GLWECtParams* params, 
                  PolyBivDFT* res, int64_t res_sl,
                  PolyBivDFT* a, int64_t a_sl,
                  PolyBivDFT* b, int64_t b_sl 
);

int normal_random_biv_poly_dft(MODULE* module, 
                               GLWECtParams*  params, 
                               PolyBivDFT* res
);


//! COMMON PART (begin)

int64_t poly_univ_bytes(GLWECtParams* params);
int64_t poly_biv_bytes(GLWECtParams* params);
void biv_to_univ(GLWECtParams* params, double* res_univ, PolyBiv* poly);


/**
 * @brief Test poly_biv_size
 */
Test(poly_biv_bytes, basic){
    GLWECtParams* params = malloc(sizeof(GLWECtParams));
    new_glwe_ct_params(params, 4, 1, 2, 8);  // N=4, k=1, kappa=2, n_limbs=8

    cr_assert_eq(poly_biv_size(params), 4, "poly_biv_size failed: got %lld, expected 4", (long long)poly_biv_size(params));
    
    delete_glwe_ct_params(params);
}

/**
 * @brief Test poly_biv_size
 */
Test(poly_biv_coef, basic){
    GLWECtParams* params = malloc(sizeof(GLWECtParams));
    new_glwe_ct_params(params, 4, 1, 2, 8);  // N=4, k=1, kappa=2, n_limbs=8

    cr_assert_eq(poly_biv_size(params), 4, "poly_biv_size failed: got %lld, expected 4", (long long)poly_biv_size(params));
    
    delete_glwe_ct_params(params);
}

/**
 * @brief Test poly_biv_bytes and poly_univ_bytes
 */
Test(poly_bytes, basic){
    GLWECtParams* params = malloc(sizeof(GLWECtParams));
    new_glwe_ct_params(params, 5, 1, 1, 10); // N=5, k=1, kappa=1, n_limbs=10

    cr_assert_eq(poly_biv_bytes(params), poly_biv_coef_number(params) * sizeof(int64_t),
                 "poly_biv_bytes failed: got %lld, expected %lld",
                 (long long)poly_biv_bytes(params), (long long)(poly_biv_coef_number(params) * sizeof(int64_t)));

    cr_assert_eq(poly_univ_bytes(params), params->N * sizeof(int64_t),
                 "poly_univ_bytes failed: got %lld, expected %lld",
                 (long long)poly_univ_bytes(params), (long long)(params->N * sizeof(int64_t)));

    delete_glwe_ct_params(params);
}

/**
 * @brief Test add_biv_poly_dft correctness with random normal DFT polynomials
 */
Test(add_biv_poly_dft, add_biv_poly_dft){
    GLWECtParams* params = malloc(sizeof(GLWECtParams));
    new_glwe_ct_params(params, 2, 1, 1, 4);
    int64_t size = poly_biv_size(params);
    int64_t total = size * params->N;

    PolyBivDFT* a = malloc(sizeof(PolyBivDFT) * total);
    PolyBivDFT* b = malloc(sizeof(PolyBivDFT) * total);
    PolyBivDFT* res = malloc(sizeof(PolyBivDFT) * total);

    MODULE* module = new_module_info(params->N, FFT64);

    normal_random_biv_poly_dft(module, params, a);
    normal_random_biv_poly_dft(module, params, b);

    add_biv_poly_dft(params, res, params->N, a, params->N, b, params->N);

    for(int64_t i = 0; i < size; i++){
        for(int64_t j = 0; j < params->N; j++){
            int64_t idx = j + i * params->N;
            cr_assert_float_eq(res[idx], a[idx] + b[idx], 1e-9,
                "add_biv_poly_dft mismatch at index %lld: %f + %f = %f, got %f",
                (long long)idx, a[idx], b[idx], a[idx] + b[idx], res[idx]);
        }
    }

    free(a); free(b); free(res);
    delete_module_info(module);
    delete_glwe_ct_params(params);
}