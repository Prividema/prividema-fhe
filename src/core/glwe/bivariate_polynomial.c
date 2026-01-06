#include "bivariate_polynomial.h"

//! BIV POLY PART (begin) 

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
 * @brief The number of coefficient in bivariate polynomial. 
 * 
 * @param params The GLWE parameters. 
 * @return int64_t 
 * 
 * @note The number of coefficient is the same in and out DFT space. 
 */
int64_t poly_biv_coef_number(GLWECtParams* params){
    int64_t N = params->N;
    return poly_biv_size(params) * N;
}

/**
 * @brief Adds two bivariate polynomial and puts it in res. 
 * 
 * @param params The GLWE parameters.
 * @param res The result bivariate polynomial.
 * @param res_sl The stride between each ZnX polynomial.
 * @param a The left-hand side bivariate polynomial.
 * @param a_sl The stride between each ZnX polynomial.
 * @param b The right-hand side bivariate polynomial.
 * @param b_sl The stride between each ZnX polynomial.
 */
void add_biv_poly(GLWECtParams* params, 
                  PolyBiv* res, int64_t res_sl,
                  PolyBiv* a, int64_t a_sl,
                  PolyBiv* b, int64_t b_sl 
){
    for(int64_t i = 0 ; i < poly_biv_size(params) ; i++)
    {
        for(int64_t p = 0 ; p < params->N ; p++)
        {
        res[p + i*res_sl] = a[p + i*a_sl] + b[p + i*b_sl];
        }
    }
    
}


//! BIV POLY IN DFT PART (begin) 

/**
 * @brief The number of coefficient in bivariate polynomial. 
 * 
 * @param params The GLWE parameters. 
 * @return int64_t 
 * 
 * @note The number of independent coefficients of a polynomial in DFT space is half the number of coefficients in ZnX, 
 * due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
int64_t poly_biv_coef_number_dft(GLWECtParams* params){
    int64_t N = params->N;
    return (poly_biv_size(params) * N)/2;
}

/**
 * @brief Adds two bivariate polynomial and puts it in res in DFT space. 
 * 
 * @param params The GLWE parameters.
 * @param res The result bivariate polynomial in DFT space.
 * @param res_sl The stride between each ZnX polynomial.
 * @param a The left-hand side bivariate polynomial in DFT space.
 * @param a_sl The stride between each ZnX polynomial.
 * @param b The right-hand side bivariate polynomial in DFT space.
 * @param b_sl The stride between each ZnX polynomial.
 */
void add_biv_poly_dft(GLWECtParams* params, 
                      PolyBivDFT* res, int64_t res_sl,
                      PolyBivDFT* a, int64_t a_sl,
                      PolyBivDFT* b, int64_t b_sl 
){
    for(int64_t i = 0 ; i < poly_biv_size(params) ; i++)
    {
        for(int64_t p = 0 ; p < params->N ; p++)
        {
        res[p + i*res_sl] = a[p + i*a_sl] + b[p + i*b_sl];
        }
    }
    
}


//! COMMON PART (begin)

/**
 * @brief The number of bytes needed to store a bivariate polynomial.
 * 
 * @param params The GLWE parameters.
 * @return int64_t 
 * 
 * @note The number of bytes needed to store a bivariate polynomial is the same in and out of DFT space.
 */
int64_t poly_biv_bytes(GLWECtParams* params){
    return poly_biv_coef_number(params) * sizeof(int64_t);
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