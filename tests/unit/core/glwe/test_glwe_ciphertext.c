#include "core/glwe/glwe_ciphertext.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "vec_znx_arithmetic_private.h"
#include "spqlios_alias.h"
#include "distributions.h"
#include "utils.h"

#define NBASE 4
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE 1e-7


//! COMMON PART (begin)

/**
 * @brief Tests whether glwe_size computes the right size of a GLWE ciphertext.
 */
Test(glwe_size, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    
    cr_assert(eq(i64, glwe_size(params), NLIMBSBASE));

    delete_glwe_ct_params(params);
}

/**
 * @brief Tests whether glwe_bytes computes the right number of bytes in a GLWE ciphertext.
 */
Test(glwe_bytes, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    
    cr_assert(eq(i64, glwe_bytes(params), NLIMBSBASE*NBASE*8));

    delete_glwe_ct_params(params);
}

/**
 * @brief Tests whether mult_vec_znx_dft multiply correctly two Zn[X] polynomials a and b. Ie res = a*b.
 */
Test(mult_vec_znx_dft, size_equal_one){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    int64_t* res = calloc(poly_univ_bytes(params),1);
    
    // def a = 1 + X
    int64_t* a = calloc(poly_univ_bytes(params), 1);
    new_uniform_random_vec(NBASE, a, 1, NBASE, 14);

    //def b = 1 + X
    int64_t* b = calloc(poly_univ_bytes(params), 1);
    new_uniform_random_vec(NBASE, b, 1, NBASE, 14);

    double* res_dft = calloc(poly_univ_bytes(params),1);
    double* a_dft = calloc(poly_univ_bytes(params),1);
    double* b_dft = calloc(poly_univ_bytes(params),1);

    vec_znx_dft_p(module, res_dft, 1, res, 1, NBASE);
    vec_znx_dft_p(module, a_dft, 1, a, 1, NBASE);
    vec_znx_dft_p(module, b_dft, 1, b, 1, NBASE);
    
    // res_dft = DFT(a*b)
    mult_vec_znx_dft(module, res_dft, 1, a_dft, 1, b_dft, 1);
    
    // res = a*b
    vec_znx_idft_p(module, res, 1, res_dft, 1);

    // Compare the real coefficient res_p for each p in [0, NBASE -1] with the res_p mult_vec_znx_dft computed coefficient.
    for(int64_t p = 0 ; p < NBASE ; p++)
    {
        int64_t acc = 0;
        for(uint64_t k = 0 ; k <= p; k++)
        {
            acc += a[k]*b[p-k];
        } 
        for(uint64_t k = p + 1; k < NBASE; k++)
        {
            acc += -a[k]*b[NBASE + p - k];
        }  
        cr_log_info("acc %ld res[p] %ld X^%ld", acc, res[p], p);
        cr_assert(eq(i64, res[p], acc));
    }

    free(res); free(res_dft); free(a); free(a_dft); free(b); free(b_dft);
    delete_module_info(module);
    delete_glwe_ct_params(params);
}

/**
 * @brief Tests whether mult_vec_znx_dft correctly multiplies two Zn[X] vectors a and b component-wise..
 * It draws a random uniform size, ie a random uniform number of Zn[X] polynomials. 
 * Ie a = (a_i), b = (b_i) -> res = (a_i * b_i). Where a_i and b_i are in Zn[X]
 */
Test(mult_vec_znx_dft, random_size){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    int64_t size = 0;

    while(size <= 0){
        rand_uniform(&size, 8);
    }

    int64_t* res = calloc(poly_univ_bytes(params)*size,1);
    
    // def a = 1 + X
    int64_t* a = calloc(poly_univ_bytes(params)*size, 1);
    new_uniform_random_vec(NBASE, a, size, NBASE, 14);

    //def b = 1 + X
    int64_t* b = calloc(poly_univ_bytes(params)*size, 1);
    new_uniform_random_vec(NBASE, b, size, NBASE, 14);

    double* res_dft = calloc(poly_univ_bytes(params)*size,1);
    double* a_dft = calloc(poly_univ_bytes(params)*size,1);
    double* b_dft = calloc(poly_univ_bytes(params)*size,1);

    vec_znx_dft_p(module, res_dft, size, res, size, NBASE);
    vec_znx_dft_p(module, a_dft, size, a, size, NBASE);
    vec_znx_dft_p(module, b_dft, size, b, size, NBASE);
    
    // res_dft = DFT(a*b)
    mult_vec_znx_dft(module, res_dft, size, a_dft, size, b_dft, size);
    
    // res = a*b
    vec_znx_idft_p(module, res, size, res_dft, size);

    // Compare the real coefficient res_p for each p in [0, NBASE -1] with the res_p mult_vec_znx_dft computed coefficient.
    for(int64_t i = 0 ; i < size ; i++){
        for(int64_t p = 0 ; p < NBASE ; p++)
        {
            int64_t acc = 0;
            for(uint64_t k = 0 ; k <= p; k++)
            {
                acc += a[i*NBASE + k]*b[i*NBASE + p-k];
            } 
            for(uint64_t k = p + 1; k < NBASE; k++)
            {
                acc += -a[i*NBASE + k]*b[i*NBASE + NBASE + p - k];
            }  
            cr_log_info("acc %ld res[i*NBASE + p] %ld X^%ld", acc, res[i*NBASE + p], p);
            cr_assert(eq(i64, res[i*NBASE + p], acc));
        }
    }

    free(res); free(res_dft); free(a); free(a_dft); free(b); free(b_dft);
    delete_module_info(module);
    delete_glwe_ct_params(params);
}

//! GLWE PART (begin)

/**
 * @brief Tests whether glwe_bytes computes the right number of coefficient in a GLWE ciphertext.
 */
Test(glwe_coef_number, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    
    cr_assert(eq(i64, glwe_coef_number(params), NLIMBSBASE*NBASE));

    delete_glwe_ct_params(params);
}

/**
 * @brief Tests whether new_glwe returns a non-NULL pointer.
 */
Test(new_glwe, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GLWECiphertext* ct = new_glwe(params);
    
    cr_assert(eq(int, (ct != NULL)&&(ct->vec != NULL), 1));

    delete_glwe(ct);
    delete_glwe_ct_params(params);
}

//! GLWE IN DFT PART (begin)

/**
 * @brief Tests whether glwe_coef_number_dft computes the right number of coefficient in a GLWE ciphertext in DFT space.
 */
Test(glwe_coef_number_dft, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

    cr_assert(eq(i64, glwe_coef_number_dft(params), NLIMBSBASE*NBASE/2));

    delete_glwe_ct_params(params);
}

/**
 * @brief Tests whether new_glwe_dft returns a non-NULL pointer.
 */
Test(new_glwe_dft, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GLWEPreparedCt* ct = new_glwe_dft(params);
    
    cr_assert(eq(int, (ct != NULL)&&(ct->pvec != NULL), 1));

    delete_glwe_dft(ct);
    delete_glwe_ct_params(params);
}


