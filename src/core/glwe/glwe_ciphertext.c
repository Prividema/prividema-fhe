#include "glwe_ciphertext.h"
#include "spqlios_alias.h"
#include "vec_znx_arithmetic_private.h"
#include <string.h>


//! GLWE PART (begin)


//! GLWE IN DFT PART (begin)


//! COMMON PART (begin)

/**
 * @brief Return the size of a bivGLWE ciphertext, in DFT space & out of DFT space.
 * 
 * @param params The GLWE parameters.
 * @return int64_t 
 * 
 * @note The size of a bivGLWE ciphertext is the same in and out of DFT space.
 */
int64_t glwe_size(GLWECtParams* params){
    return params->n_limbs;
}

/**
 * @brief The number of bytes needed to store a bivGLWE ciphertext.
 * 
 * @param params The GLWE parameters.
 * @return int64_t 
 * 
 * @note The number of bytes needed to store a bivGLWE ciphertext, is the same in and out of DFT space. 
 */
int64_t glwe_bytes(GLWECtParams* params){
    int64_t N = params->N;
    return glwe_size(params) * N * sizeof(int64_t); 
}

/**
 * @brief The size of a bivariate polynomial.
 * 
 * @param params 
 * @return int64_t 
 * 
 * @note The size of a bivariate polynomial is the same in and out of DFT space.
 */
int64_t poly_biv_size(GLWECtParams* params){
    return params->n_limbs/(params->k + 1);
}

/**
 * @brief The number of bytes needed to store a bivariate polynomial.
 * 
 * @param params The GLWE parameters.
 * @return int64_t 
 * 
 * @note The number of bytes needed to store a bivariate polynomial is the same in and out of DFT space.
 */
int64_t poly_biv_bytes(GLWECtParams* params){
    int64_t N = params->N;
    return poly_biv_size(params) * N * sizeof(int64_t);
}

/**
 * @brief The number of bytes needed to store a univariate polynomial.
 * 
 * @param params The GLWE parameters.
 * @return int64_t 
 * 
 * @note The number of bytes needed to store a univariate polynomial is the same in and out of DFT space.
 */
int64_t poly_univ_bytes(GLWECtParams* params){
    int64_t N = params->N;
    return N * sizeof(int64_t);
}

/**
 * @brief Compute the polynomial product of c and d, component-wise in DFT space.
 * 
 * @param module The module stocking the degree N.
 * @param res_dft The result in DFT space.
 * @param res_size The result's size.
 * @param c_dft The left-hand side polynomial in DFT space .
 * @param c_size The left-hand size of c_dft.
 * @param d_dft The right-hand side polynomial in DFT space.
 * @param d_size The right-hand size of c_dft.
 * 
 * @note `res_dft = ( DFT(c_0) * DFT(d_0) , ... , DFT(c_smin) * DFT(d_smin) , 0's)`. There are enough 0's to match the size of res_dft.
 */
void vec_znx_dft_mult(const MODULE* module, 
              double* res_dft, int64_t res_size,
              double* c_dft, int64_t c_size,  
              double* d_dft, int64_t d_size
){
    int64_t N = module->nn;

    if (c_size < d_size){
        int64_t smin = c_size < res_size ? c_size : res_size;
        
        for (int i = 0 ; i < smin; i++){
            for (int64_t j = 0 ; j < N/2 ; j++){ 
                // i*N + j corresponds to the p-th coefficient's Im[DFT(c)] & Im[DFT(d)]
                double c_re = ((double*)c_dft)[i*N + j];
                double d_re = ((double*)d_dft)[i*N + j];

                // i*N + j + N/2 corresponds to the p-th coefficient's Im[DFT(c)] & Im[DFT(d)]
                double c_im = ((double *)c_dft)[i*N + j + N/2];
                double d_im = ((double *)d_dft)[i*N + j];

                ((double*)res_dft)[i*N + j] = c_re * d_re - c_im * d_im;
                ((double*)res_dft)[i*N + j + N/2] = c_re * d_im + c_im * d_re; 
            }
        }
        
        // fill up remaining part with 0's
        double* const dres_dft = (double*)res_dft;
        memset(dres_dft + smin* N, 0, (res_size - smin) * N * sizeof(double));
    }
    else {
        vec_znx_dft_mult(module, res_dft, res_size, d_dft, d_size, c_dft, c_size);
    }
}

/**
 * @brief Compute the univariate repreentation in RnX of a bivariate polynomial.
 * 
 * @param params The GLWE parameters.
 * @param res_univ The result univariate polynomial.
 * @param poly The input bivariate polynomial.
 */
void biv_to_univ(GLWECtParams* params, double* res_univ, PolyBiv* poly){
    // GLWE parameters 
    int64_t N = params->N; 
    int64_t k = params->k; 
    int64_t l = params->n_limbs / (k+1);

    // acc_(i+1) = acc_i + limb_i(poly) * 1/Bg^i
    double* acc = malloc(N * sizeof(double));
    for(int64_t i = 0 ; i < l ; i++){
        for(int64_t p = 0 ; p < N ; p++){
            acc[p] = (double)poly[N*i + p] / ((double) params->kappa)
        }
    }
}