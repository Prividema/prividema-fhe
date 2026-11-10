#ifndef bivGLWE_CT_PARAMS_H
#define bivGLWE_CT_PARAMS_H

#include <stdint.h>

/**
 * @file glwe_params.h
 *
 * @brief Parameters object for a GLWE problem/ciphertext.
 *
 * This file contains the definition of GLWE parameters, as well as accessor/helper functions for parameters
 * that can be derived from them.
 */

/**
 * @brief enum for the available distributions to sample the error/noise used in FHE encryption operations
 */
typedef enum noise_type_t {
	NOISE_UNIFORM_POWER_OF_TWO,  ///< Uniform noise, with the range being restricted to powers of 2. The smallest suitable one is selected from sigma.
	NOISE_UNIFORM,   ///< Uniform noise in [-B, B] on the last limb, B chosen from sigma (no power-of-2 restriction).
	                 ///< This is slower (rejection sampling), albeit more precise
	NOISE_BINOMIAL,  ///< Centered binomial noise on the last limb (B(n, 1/2) - n/2, n chosen from sigma).
	                 ///< Costs O(n) random bits per coefficient: only for small noise
	NOISE_NORMAL     ///< Rounded normal noise on the last limb (Box-Muller, not constant-time)
} NoiseType;

/**
 * @brief Parameters object for a GLWE problem/ciphertext
 */
typedef struct glwe_ct_params
{
	uint64_t nn;                   ///< The cyclotomic polynomial \cyclopol degree
	uint64_t k;                    ///< Number of secret key terms, \f$ k=1 \f$ for RLWE
	uint64_t kappa;                ///< The \K for the base-\f$2^K\f$ representation
	uint64_t ciphertext_nb_limbs;  ///< @brief Number of limbs across all polynomials.
	                               ///<
	                               ///< In other words, \f$ \ell_a \cdot k + \ell_b \f$.
	                               ///<
	///< Only \f$ \ell_a = \ell_b \f$ and \f$ \ell_a = \ell_b + 1 \f$ are supported due to memory layout.
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
 * @param n_limbs \f$ \ell_a \cdot k + \ell_b \f$
 * @param sigma The standard deviation of the error distribution, on the torus scale (e.g. 2^-142).
 *              It is internally rescaled by \f$ 2^{K \cdot \ell_a} \f$ and must then fit in [1, 2^63].
 *              A value of 0 disables the noise.
 * @param noise_type The type of noise that will be used for encryption
 * @return The newly allocated parameters object, or NULL if the allocation failed, if the rescaled sigma
 *         does not fit (too big or too small), if the noise type is not supported yet, or if the parameters are
 *         invalid: nn == 0, k == 0, kappa outside [1, 63], sigma negative or not finite, or n_limbs not of the
 *         form k * l_a + l_b with l_a >= 1 and l_b in {l_a - 1, l_a}
 *
 */
GLWEParams* new_glwe_params(uint64_t nn, uint64_t k, uint64_t kappa, uint64_t n_limbs, double sigma,
                            NoiseType noise_type);

/**
 * @brief Deletes the set of bivGLWE parameters.
 *
 * @param params_glwe The GLWE parameters.
 */
void delete_glwe_params(GLWEParams* params_glwe);

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
 * In other words, \f$ \ell_a \cdot k + \ell_b \f$
 *
 * @param params_glwe The GLWE parameters.
 */
uint64_t glwe_params_n_limbs(const GLWEParams* params_glwe);

/**
 * @brief Returns the standard deviation of the noise that is added for encryption
 *
 * @param params_glwe The GLWE parameters.
 *
 * @return The standard deviation of the noise, on the torus scale:
 *         - NOISE_NORMAL: sigma as given to new_glwe_params;
 *         - NOISE_UNIFORM: sqrt(B(B+1)/3) * 2^(-K l_a), B being uniform_range;
 *         - NOISE_BINOMIAL: sqrt(binomial_n)/2 * 2^(-K l_a);
 *         - NOISE_UNIFORM_POWER_OF_TWO: NAN (the power of 2 only bounds the noise).
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
 * @brief The number of coefficients in a bivariate GLWE ciphertext in the DFT domain.
 *
 * @param params_glwe The bivariate GLWE parameters.
 * @return The number of coefficients in a bivariate GLWE ciphertext in the DFT domain.
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * \ZnX, due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t glwe_coef_number_dft(const GLWEParams* params_glwe);

/**
 *
 * @brief Returns the smallest difference representable in the bivariate
 * decomposition described by the provided params, i.e. 2^(-K * l_b)
 *
 * @param params_glwe The GLWE parameters.
 */
double glwe_bivariate_epsilon(const GLWEParams* params_glwe);

#endif  // bivGLWE_CT_PARAMS_H
