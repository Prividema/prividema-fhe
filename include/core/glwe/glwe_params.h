#ifndef bivGLWE_CT_PARAMS_H
#define bivGLWE_CT_PARAMS_H

#include <stdint.h>

typedef struct glwe_ct_params
{
	uint64_t nn;                   // The cyclotomic polynomial degree
	uint64_t k;                    // Number of secret key terms, k=1 for RLWE
	uint64_t kappa;                // Used for the base-2^kappa representation (K)
	uint64_t ciphertext_nb_limbs;  // Number of limbs across all polynomials. In other words, l_a*k + l_b
	double sigma;                  // The standard deviation of the error distribution.
} GLWEParams;

/**
 * @brief Creates a new set of bivGLWE parameters.
 *
 * @param nn The polynomial maximum degree in X
 * @param k The number of Zn[X] polynomial in the secret key.
 * @param kappa The exponent for the base-2^kappa representation.
 * @param n_limbs l_a * k + l_b
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
 * @brief Returns the parameter l_a (l for the As in the GLWE)
 *
 * @param params_glwe
 * @return uint64_t
 *
 */
uint64_t glwe_params_l_a(const GLWEParams* params_glwe);

/**
 * @brief Returns the parameter l_b (l for the b element of the GLWE)
 *
 * @param params_glwe
 * @return uint64_t
 *
 */
uint64_t glwe_params_l_b(const GLWEParams* params_glwe);

/**
 * @brief Returns the number of limbs across all polynomials in a GLWE ciphertext
 * In other words, l_a * k + l_b
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

double glwe_bivariate_epsilon(const GLWEParams* params);

#endif  // bivGLWE_CT_PARAMS_H
