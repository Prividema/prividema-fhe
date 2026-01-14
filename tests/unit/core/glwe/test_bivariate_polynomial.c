#include "bivariate_polynomial.h"
#include "core/glwe/glwe.h"
#include "common/distributions.h"
#include "vec_znx_arithmetic_private.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <math.h>

#define NBASE 2
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)

//! COMMON PART (begin)

/**
 * @brief Test poly_biv_size
 */
Test(poly_biv_size, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);

    cr_assert(eq(i64, poly_biv_size(params), LBASE, 
                "poly_biv_size failed: got %lld, expected %lld", poly_biv_size(params), LBASE));
    
    delete_glwe_ct_params(params);
}

/**
 * @brief Test poly_biv_bytes
 */
Test(poly_univ_bytes, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);

    cr_assert(eq(i64, poly_univ_bytes(params), NBASE*sizeof(int64_t), 
                "poly_univ_bytes failed: got %lld, expected %lld", poly_univ_bytes(params), NBASE*sizeof(int64_t)));
    
    delete_glwe_ct_params(params);
}

/**
 * @brief Test poly_biv_bytes
 */
Test(poly_biv_bytes, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);

    cr_assert(eq(i64, poly_biv_bytes(params), NBASE*LBASE*sizeof(int64_t), 
                "poly_biv_bytes failed: got %lld, expected %lld", poly_biv_bytes(params), NBASE*LBASE*sizeof(int64_t)));
    
    delete_glwe_ct_params(params);
}

/**
 * @brief Test biv_to_univ
*/
Test(biv_to_univ, test_with_random_biv_generation ){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    double* pol_univ = malloc(poly_univ_bytes(params));
    PolyBiv* pol_biv = new_normal_random_biv_poly(module, params);
    
    biv_to_univ(params, pol_univ, pol_biv);

    // for(int64_t p = 0 ; p < NBASE ; p++){
    //     cr_log_info("%e X^%ld", pol_univ[p], p);
    // }

    free(pol_univ); free(pol_biv);
    delete_glwe_ct_params(params);
    delete_module_info(module);
} 

/**
 * @brief Test univ_to_biv
 */
Test(univ_to_biv, one_test){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    double* pol_univ = malloc(poly_univ_bytes(params));
    normal_random_vec(NBASE, pol_univ, 1, NBASE, 0.0, 1e-2);

    int64_t mask = (1LL << KAPPABASE) - 1;
    for(int64_t p = 0 ; p < NBASE ; p++){
        cr_log_info("A %e X^%ld", pol_univ[p], p);
        for(int64_t i = 1 ; i < LBASE ; i++){
            // cr_log_info("A(XY) %e Y^%ld", ldexp(pol_univ[p], i*KAPPABASE), i) ;
            // cr_log_info("A(XY) %ld Y^%ld", (int64_t) ldexp(pol_univ[p], i*KAPPABASE) & mask, i) ;
            }
    }

    PolyBiv* pol_biv = malloc(poly_biv_bytes(params));
    univ_to_biv(params, pol_biv, pol_univ);
    
    for(int64_t p = 0 ; p < NBASE ; p++){
        cr_log_info("A %e X^%ld", pol_univ[p], p);
        for(int64_t i = 1 ; i < LBASE ; i++)
            cr_log_info("A(XY) %ld Y^%ld", pol_biv[i*NBASE + p], i);
    }

    for(int64_t p = 0 ; p < NBASE ; p++){
        double acc = 0; 
        for(int64_t i = 1 ; i < LBASE ; i++){
            acc += ldexp((double)pol_biv[i * NBASE + p], -i * KAPPABASE);
        }
        cr_log_info("acc %lf pol %lf p %ld", acc, pol_univ[p], p);
        cr_assert(epsilon_eq(dbl, acc - pol_univ[p], 0, ldexp(1.0,-(LBASE-1)*KAPPABASE)), "acc %e pol %e p %ld", acc, pol_univ[p], p);
        
    }
    free(pol_univ); free(pol_biv);
    delete_glwe_ct_params(params);
    delete_module_info(module);
}

/**
 * @brief Test univ_to_biv
 
Test(univ_to_biv, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    double* pol_univ = malloc(poly_univ_bytes(params));
    normal_random_vec(NBASE, pol_univ, 1, NBASE, 0.0, 1e-2);

    int64_t mask = (1LL << KAPPABASE) - 1;
    for(int64_t p = 0 ; p < NBASE ; p++){
        cr_log_info("A %e X^%ld", pol_univ[p], p);
        for(int64_t i = 1 ; i < LBASE ; i++){
            // cr_log_info("A(XY) %e Y^%ld", ldexp(pol_univ[p], i*KAPPABASE), i) ;
            // cr_log_info("A(XY) %ld Y^%ld", (int64_t) ldexp(pol_univ[p], i*KAPPABASE) & mask, i) ;
            }
    }

    PolyBiv* pol_biv = malloc(poly_biv_bytes(params));
    univ_to_biv(params, pol_biv, pol_univ);
    
    for(int64_t p = 0 ; p < NBASE ; p++){
        cr_log_info("A %e X^%ld", pol_univ[p], p);
        for(int64_t i = 1 ; i < LBASE ; i++)
            cr_log_info("A(XY) %ld Y^%ld", pol_biv[i*NBASE + p], i);
    }

    for(int64_t p = 0 ; p < NBASE ; p++){
        double acc = 0; 
        for(int64_t i = 1 ; i < LBASE ; i++){
            acc += ldexp((double)pol_biv[i * NBASE + p], -i * KAPPABASE);
        }
        cr_log_info("acc %lf pol %lf p %ld", acc, pol_univ[p], p);
        cr_assert(epsilon_eq(dbl, acc - pol_univ[p], 0, ldexp(1.0,-(LBASE-1)*KAPPABASE)), "acc %e pol %e p %ld", acc, pol_univ[p], p);
        
    }
    free(pol_univ); free(pol_biv);
    delete_glwe_ct_params(params);
    delete_module_info(module);
}*/

/**
 * @brief Test univ_to_biv
 
Test(univ_to_biv, maths_test){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    double* pol_univ = malloc(poly_univ_bytes(params));
    normal_random_vec(NBASE, pol_univ, 1, NBASE, 0.0, 1e-7);

    PolyBiv* pol_biv = malloc(poly_biv_bytes(params));
    univ_to_biv(params, pol_biv, pol_univ);
    
    double* pol_univ_bis = malloc(poly_univ_bytes(params));
    biv_to_univ(params, pol_univ_bis, pol_biv);

    for(int64_t p = 0 ; p < NBASE ; p++){
        cr_log_info("A %e X^p", pol_univ[p]);
        cr_log_info("A' %e X^p", pol_univ_bis[p]);
    }
    
    free(pol_univ); free(pol_univ_bis); free(pol_biv);
    delete_glwe_ct_params(params);
    delete_module_info(module);
}*/


//! BIV POLY PART (begin) 

/**
 * @brief Tests if it returns the right size for N in [1,100]
 * @note n_limbs = (k + 1) * l
 */
Test(poly_biv_coef_number, classic_params){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE); 

    cr_assert(eq(i64, poly_biv_coef_number(params), NBASE * LBASE, 
                "poly_biv_coef_number failed: got %lld, expected %lld", poly_biv_coef_number(params), NBASE * LBASE));

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
 * @brief Test new_normal_random_biv_poly
 * 
 */
Test(new_normal_random_biv_poly, is_it_working){
    MODULE* module = new_module_info(NBASE, FFT64);
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    PolyBiv* a = new_normal_random_biv_poly(module, params);

    for(int64_t i = 0 ; i < LBASE ; i++)
    {
        for(int64_t p = 0 ; p < params->N ; p++)
        {
            cr_assert(le(i64, a[i*NBASE + p], (1 << (KAPPABASE-1)), "The coefficient of a(X^p, Y^i) is greater than 2^(kappa-1)."));
            cr_assert(ge(i64, a[i*NBASE + p], -(1 << (KAPPABASE-1)), "The coefficient of a(X^p, Y^i) is smaller than -2^(kappa-1)."));
        }
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

    for(int64_t i = 0; i < LBASE; i++)
    {
        for(int64_t p = 0; p < NBASE; p++)
        {
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

    
    for(int64_t i = 0; i < LBASE; i++)
    {
        for(int64_t p = 0; p < NBASE; p++)
        {
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

/**
 * @brief Test new_normal_random_biv_poly_dft
 * 
 */
Test(new_normal_random_biv_poly_dft, is_it_working){
    MODULE* module = new_module_info(NBASE, FFT64);
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    
    PolyBivDFT* a_dft = new_normal_random_biv_poly_dft(module, params);
    PolyBiv* a = malloc(poly_biv_bytes(params));
    PolyBiv* a_normalized = malloc(poly_biv_bytes(params));

    vec_znx_idft(module, (VEC_ZNX_BIG *)a, LBASE, (VEC_ZNX_DFT *)a_dft, LBASE, NULL);
    uint8_t* tmp = malloc(vec_znx_normalize_base2k_tmp_bytes(module));
    vec_znx_normalize_base2k(module, KAPPABASE, a_normalized, LBASE, NBASE, a, LBASE, NBASE, tmp);

    for(int64_t i = 0 ; i < LBASE ; i++)
    {
        for(int64_t p = 0 ; p < params->N ; p++)
        {
            cr_assert(le(i64, a_normalized[i*NBASE + p], (1 << (KAPPABASE-1)), "The coefficient of a(X^p, Y^i) is greater than 2^(kappa-1)."));
            cr_assert(ge(i64, a_normalized[i*NBASE + p], -(1 << (KAPPABASE-1)), "The coefficient of a(X^p, Y^i) is smaller than -2^(kappa-1)."));
        }
    }
    cr_assert(1);

    free(tmp);
    free(a_dft);
    free(a);
    free(a_normalized);
    delete_glwe_ct_params(params);
    delete_module_info(module);
}