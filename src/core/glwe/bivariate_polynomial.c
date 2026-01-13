#include "bivariate_polynomial.h"
#include "distributions.h"
#include <math.h>
#include <stdio.h>
#include <string.h>


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
                                    GLWECtParams* params
){
    // GLWE parameters
    int64_t N = params->N;
    int64_t kappa = params->kappa;
    int64_t l = params->n_limbs/(params->k + 1);

    // Base-2Kappa normalized bivariate polynomial in ZnXY
    PolyBiv* pol = malloc(poly_biv_bytes(params));
    if(pol==NULL){
        perror("Malloc failed.");
        return NULL;
    }
    // Univariate polynomial in RnX
    double* tmp_pol_inR_univ = malloc(poly_univ_bytes(params));
    if(tmp_pol_inR_univ==NULL){
        perror("Malloc failed.");
        free(pol);
        return NULL;}

    // Bivariate polynomial in ZnXY
    int64_t* tmp_pol_inZ = malloc(poly_biv_bytes(params));
    if(tmp_pol_inZ==NULL){
        perror("Malloc failed.");
        free(pol);
        free(tmp_pol_inR_univ);
        return NULL;}
    
    // Fills each tmp_pol_inZ(X^p, Y^i) with coefficients in [-2^(kappa* - 1) ; 2^(kappa - 1) - 1]
    int64_t mask = (1 << (kappa + 1)) - 1;
    for(int64_t p = 0 ; p < N ; p++){
        if(rand_normal(tmp_pol_inR_univ + p, 0.0, 1e-7) < 0) 
                return NULL;
                // TODO Test
        for(int64_t i = 0 ; i < l ; i++){
            // tmp_pol_inZ(X^p, Y^i) = the i_ème block of kappa bits, starting from the MSB, of tmp_pol_inR_univ(X^p)
            tmp_pol_inZ[i*N + p] = ((int64_t)ldexp(tmp_pol_inR_univ[p], i*kappa)) & mask;
        }
    }

    // Then does a base-2Kappa normalization 
    vec_znx_normalize_base2k_p(module, kappa, pol, poly_biv_size(params), N, tmp_pol_inZ, poly_biv_size(params), N);

    free(tmp_pol_inR_univ);
    free(tmp_pol_inZ);

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
 * @return int 
 */
PolyBivDFT* new_normal_random_biv_poly_dft(MODULE* module, 
                                           GLWECtParams* params
){
    // Base-2Kappa normalized bivariate polynomial in DFt space
    PolyBivDFT* pol_dft = malloc(poly_biv_bytes(params));
    if(pol_dft == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    // Base-2Kappa normalized bivariate polynomial
    PolyBiv* tmp_pol = new_normal_random_biv_poly(module, params);
    if(tmp_pol == NULL){
        perror("Malloc failed.");
        free(pol_dft);
        return NULL;
    }
    
    // Then compute in DFT space
    vec_znx_dft_p(module, pol_dft, poly_biv_size(params), tmp_pol, poly_biv_size(params), params->N);

    free(tmp_pol);
    
    return pol_dft;
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

/**
 * @brief Computes P(X,2^(-kappa)) for P a bivariate polynomial.
 * 
 * @param params The GLWE parameters.
 * @param res_univ The result univariate polynomial in RnX.
 * @param poly The input bivariate polynomial.
 */
void biv_to_univ(GLWECtParams* params, double* res_univ, PolyBiv* poly){
    //GLWE parameters
    int64_t N = params->N;
    int64_t kappa = params->kappa;
    int64_t l = poly_biv_size(params);

    // res_univ(X^p) = ∑_i{1,l}[poly(X^p, Y^i) * 2^(-kappa*i)]
    for(int64_t i = 1 ; i < l ; i++){
        for(int64_t p = 0 ; p < N ; p++){
            res_univ[p] += ldexp((double)poly[i*N + p], - i*kappa);
        }
    }
}