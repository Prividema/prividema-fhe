#ifndef bivGLWE_CT_PARAMS_H
#define bivGLWE_CT_PARAMS_H

#include <stdint.h>

/**
 * @brief Parameters object for a GLWE problem/ciphertext
 *
 *
 */
typedef struct glwe_ct_params
{
	uint64_t nn;                   ///< The cyclotomic polynomial \cyclopol degree
	uint64_t k;                    ///< Number of secret key terms, \f$ k=1 \f$ for RLWE
	uint64_t kappa;                ///< The \K for the base-\f$2^K\f$ representation
	uint64_t ciphertext_nb_limbs;  ///< @brief Number of limbs across all polynomials.
	                               ///<
	                               ///< In other words, \f$ l_a \cdot k + l_b \f$.
	                               ///<
	                               ///< Only \l_a = \l_b and \l_a = \l_b + 1 are supported due to memory layout.
	double sigma;                  ///< The standard deviation of the error distribution.
} GLWEParams;

/**
 * @brief Creates a new set of bivariate GLWE parameters.
 *
 * @param nn The degree of the cyclotomic polynomial used for the modulo
 * @param k The number of \ZnX polynomials in a secret key.
 * @param kappa The exponent for the base-2^kappa representation.
 * @param n_limbs \f$ l_a \cdot k + l_b \f$
 * @param sigma The standard deviation of the error distribution.
 * @return The newly allocated parameters object, or NULL if it failed the allocation
 *
 */
GLWEParams* new_glwe_params(uint64_t nn, uint64_t k, uint64_t kappa, uint64_t n_limbs, double sigma);

/**
 * @brief Deletes the set of bivGLWE parameters.
 *
 * @param params The GLWE parameters.
 */
void delete_glwe_params(GLWEParams* params);

/**
 * @brief Returns the parameter \l_a (\l for the As in the GLWE)
 *
 * @param params_glwe The GLWE parameters.
 */
uint64_t glwe_params_l_a(const GLWEParams* params_glwe);

/**
 * @brief Returns the parameter \l_b (\l for the b element of the GLWE)
 *
 * @param params_glwe The GLWE parameters.
 */
uint64_t glwe_params_l_b(const GLWEParams* params_glwe);

/**
 * @brief Returns the number of limbs across all polynomials in a GLWE ciphertext
 * In other words, \f$ l_a \cdot k + l_b \f$
 *
 * @param params_glwe The GLWE parameters.
 */
uint64_t glwe_params_n_limbs(const GLWEParams* params_glwe);

/**
 * @brief The number of bytes needed to store a GLWE ciphertext in bivariate form.
 *
 * @param params_glwe The GLWE parameters.
 *
 * @note The number of bytes needed to store a bivGLWE ciphertext is the same in and out of the DFT domain.
 */
uint64_t glwe_params_bytes(const GLWEParams* params_glwe);

/**
 * @brief Return the number of coefficients in a bivariate GLWE ciphertext.
 *
 * @param params_glwe The GLWE parameters.
 */
uint64_t glwe_coef_number(const GLWEParams* params_glwe);

/**
 *
 * @brief Returns the smallest difference representable in the bivariate
 * decomposition described by the provided params
 *
 * @param params The GLWE parameters.
 */
double glwe_bivariate_epsilon(const GLWEParams* params);

#endif  // bivGLWE_CT_PARAMS_H
