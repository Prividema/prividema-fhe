#ifndef bivGLWE_CT_PARAMS_H
#define bivGLWE_CT_PARAMS_H

#include <stdint.h>

/**
 *
 * @brief enum for the available distributions to sample the error/noise used in FHE encryption operations
 *
 */
typedef enum noise_type_t {
	NOISE_UNIFORM_POWER_OF_TWO,  ///< Uniform noise, with the range being restricted to powers of 2. The smallest suitable one is selected from sigma.
	NOISE_UNIFORM,   ///< Uniform noise withouth the power-of-2 restriction in its internal sampling.
	                 ///< This will make it significantly slower, albeit potentially more precise
	NOISE_BINOMIAL,  ///< Binomial noise (on the last limb, with p=0.5, ie, symmetric)
	NOISE_NORMAL     ///< Normal noise (implementation at the moment is a proof of concept)
} NoiseType;

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
	union {
		double normal_sigma;  ///< When the noise type is normal, this contains the stdev of the normal
		uint64_t
		    binomial_n;  ///< When the noise type is binomial, this contains the parameter n used to sample a binomial on the last limb
		int64_t
		    uniform_range;  ///< When a uniform noise is selected, this value encodes the limits of the uniform distribution to be added to the last limb
		uint64_t
		    fast_uniform_nb_bits;  ///< When a "fast" uniform noise type is selected, this contains the number of bits of uniform randomness
		                           ///< that the last limb should have added
	};

	NoiseType noise_type;  ///< The type of noise that will be used for encryption
} GLWEParams;

/**
 * @brief Creates a new set of bivariate GLWE parameters.
 *
 * @param nn The degree of the cyclotomic polynomial used for the modulo
 * @param k The number of \ZnX polynomials in a secret key.
 * @param kappa The exponent for the base-2^kappa representation.
 * @param n_limbs \f$ l_a \cdot k + l_b \f$
 * @param sigma The standard deviation of the error distribution.
 * @param noise_type The type of noise that will be used for encryption
 * @return The newly allocated parameters object, or NULL if it failed the allocation
 *
 */
GLWEParams* new_glwe_params(uint64_t nn, uint64_t k, uint64_t kappa, uint64_t n_limbs, double sigma,
                            NoiseType noise_type);

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
 * @brief Returns the standard deviation of the noise that is added for encryption
 *
 * @param params_glwe The GLWE parameters.
 */
double glwe_params_stdev(const GLWEParams* params_glwe);

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
 * @brief The number of coefficient in a bivariate GLWE ciphertext in the DFT domain.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return The number of coefficient in a bivariate GLWE ciphertext in the DFT domain.
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * \ZnX, due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t glwe_coef_number_dft(const GLWEParams* params_glwe);

/**
 *
 * @brief Returns the smallest difference representable in the bivariate
 * decomposition described by the provided params
 *
 * @param params The GLWE parameters.
 */
double glwe_bivariate_epsilon(const GLWEParams* params);

#endif  // bivGLWE_CT_PARAMS_H
