#ifndef bivGLWE_CT_PARAMS_H
#define bivGLWE_CT_PARAMS_H

#include <stdint.h>
#include <stdlib.h>

typedef struct glwe_ct_params
{
	uint64_t nn;     // The cyclotomic polynomial degree, in the code m = nn / 2, ie the dimension of the DFT domain
	uint64_t k;      // Number of secret key terms, k=1 for RLWE
	uint64_t kappa;  // Used for the base-2^kappa representation (K)
	uint64_t l;      // l
	                 // at each limb we have a polynomial of degree nn (Zn[X])
	double sigma;    // The standard deviation of the error distribution.
} GLWEParams;

/**
 * @brief Creates a new set of bivGLWE parameters.
 *
 * @param nn The polynomial maximum degree in X
 * @param k The number of Zn[X] polynomial in the secret key.
 * @param kappa The exponent for the base-2^kappa representation.
 * @param n_limbs (k+1)*l.
 * @param sigma The standard deviation of the error distribution.
 * @return GLWEParams*
 */
GLWEParams* new_glwe_params(uint64_t nn, uint64_t k, uint64_t kappa, uint64_t n_limbs, double sigma);

/**
 * @brief Deletes the set of bivGLWE parameters.
 *
 * @param params The bivGLWE parameters.
 */
void delete_glwe_params(GLWEParams* params);

/**
 * @brief Returns the size - the number l of Zn[X] coefficients - for a bivariate polynomial.
 *
 * @param params_glwe
 * @return uint64_t
 *
 * @note The size of a bivariate polynomial is the same in and out of the DFT domain.
 */
uint64_t glwe_params_l(const GLWEParams* params_glwe);

/**
 * @brief Returns the number of limbs across all polynomials in a GLWE ciphertext
 * In other words, (k+1)*l
 *
 * @param params_glwe
 * @return uint64_t
 *
 */
uint64_t glwe_params_n_limbs(const GLWEParams* params_glwe);

/**
 * @brief The number of bytes needed to store a bivGLWE ciphertext.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return int64_t
 *
 * @note The number of bytes needed to store a bivGLWE ciphertext, is the same in and out of the DFT domain.
 */
uint64_t glwe_params_bytes(const GLWEParams* params_glwe);

/**
 * @brief Return the number of coefficient in a bivariate bivGLWE ciphertext.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return int64_t
 */
uint64_t glwe_coef_number(const GLWEParams* params_glwe);

#endif  // bivGLWE_CT_PARAMS_H
