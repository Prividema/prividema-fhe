#include "bivariate_polynomial.h"
#include "distributions.h"
#include <math.h>

//! BIV POLY PART (begin) 

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
 * @brief Computes a random normal bivariate polynomial.
 * 
 * @param module The module holding the degree N and FFT64.
 * @param params The GLWE parameters.
 * @param res The result bivariate polynomial. 
 * @param res_sl The result's vector stride.
 * @return int 
 */
PolyBiv* new_normal_random_biv_poly(MODULE* module, 
                                GLWECtParams*  params
){
    // GLWE parameters
    int64_t N = params->N;
    int64_t kappa = params->kappa;
    int64_t l = params->n_limbs/(params->k + 1);

    PolyBiv* pol = malloc(poly_biv_bytes(params));
    PolyBiv* tmp_pol = malloc(poly_biv_bytes(params));
    if(pol==NULL || tmp_pol == NULL)
        perror("Malloc failed.");
        return NULL;
    
    // Fills tmp_biv_pol(Y^l) with coefficients in [-2^(kappa*l - 1) ; 2^(kappa*l - 1)]
    for(int64_t p = 0 ; p < N ; p++)
    {
        if(rand_normal((double *)tmp_pol + N*l + p, 0.0, 1.0) < 0) 
                return -1;

        tmp_pol[N*l + p] = (int64_t) ldexp(tmp_pol[N*l + p], kappa * l);
    }

    // Then does a base-2Kappa normalization 
    vec_znx_normalize_base2k_p(module, kappa, pol, poly_biv_size(params), N, tmp_pol, poly_biv_size(params), N);

    return pol;
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
 * @brief Computes a random normal bivariate polynomial in DFT space.
 * 
 * @param module The module holding the degree N and FFT64.
 * @param params The GLWE parameters.
 * @param res The result bivariate polynomial in DFT space. 
 * @param res_sl The result's vector stride.
 * @return int 
 */
int normal_random_biv_poly_dft(MODULE* module, 
                               GLWECtParams*  params, 
                               PolyBivDFT* res
){
    // GLWE parameters
    int64_t N = params->N;
    int64_t kappa = params->kappa;
    int64_t l = params->n_limbs/(params->k + 1);

    int64_t* tmp_biv_pol = malloc(poly_biv_bytes(params));
    
    // Fills res(Y^l) with coefficients in [-2^(kappa*l - 1) ; 2^(kappa*l - 1)]
    for(int64_t p = 0 ; p < N ; p++)
    {
        if(rand_normal((double *)tmp_biv_pol + N*l + p, 0.0, 1.0) < 0) 
                return -1;

        tmp_biv_pol[N*l + p] = (int64_t) ldexp(tmp_biv_pol[N*l + p], kappa * l);
    }

    // Then does a base-2Kappa normalization 
    int64_t* tmp_biv_pol_normalized = malloc(poly_biv_bytes(params));

    // Then compute in DFT space
    vec_znx_dft_p(module, res, poly_biv_size(params), tmp_biv_pol_normalized, poly_biv_size(params), N);

    return 0;
}

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

