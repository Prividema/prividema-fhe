#ifndef RNG_H
#define RNG_H

#include <stddef.h>
#include <stdint.h>

#include "backend.h"
#include "maths_structures.h"

/**
 *
 * @file rng.h
 *
 * Prividema random number generation (takes randomness from system source,
 * provides functions to sample different distributions)
 */

// =============================================
// |                                           |
// |              Random Numbers               |
// |                                           |
// =============================================

/**
 * @brief Generates a uniformly sampled random number in [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param result  The resulting uniformly sampled integer
 * @param nb_bits The number of bits of the result
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int rand_uniform_pow2(int64_t* result, uint64_t nb_bits);

/**
 * @brief Generates a uniformly sampled random number in [limit_down, limit_up] via
 * power-of-2 sampling and resampling if out-of-bounds
 *
 * @param result      The resulting uniformly sampled integer
 * @param limit_down  The lower bound of the uniform sample
 * @param limit_up  The lower bound of the uniform sample
 *
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int rand_uniform(int64_t* result, int64_t limit_down, int64_t limit_up);

/**
 * @brief Generates a gaussianly sampled random number with the given parameters.
 *
 * @param result The resulting normally (Gaussian) sampled number
 * @param mu     The mean parameter for the distribution.
 * @param sigma  The standard deviation of the distribution.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int rand_normal(double* result, double mu, double sigma);

// =============================================
// |                                           |
// |      Random Polynomials and Vectors       |
// |                                           |
// =============================================

/**
 * @brief Generates a uniformly random \ZnX polynomial
 *
 * Coefficients are uniformly sampled in range [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param res     The result uniformly drawn \ZnX polynomial.
 * @param nn      Number of coeffients in the polynomial (eq. degree of the cyclotomial poly)
 * @param nb_bits Number of randomness bits per coefficient.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int uniform_pow2_random_pol_znx(PolyUniv* res, uint64_t nn, uint64_t nb_bits);

/**
 * @brief Generates a uniformly random binary \ZnX polynomial
 *
 * Coefficients are uniformly sampled in {0, 1}
 *
 * @param res     The result uniformly drawn \ZnX polynomial.
 * @param nn      Number of coeffients in the polynomial (eq. degree of the cyclotomial poly)
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int binary_random_pol_znx(PolyUniv* res, uint64_t nn);

/**
 * @brief Generates a uniformly random \ZnX polynomial
 *
 * Coefficients are uniformly sampled in range [low_bound, high_bound]
 *
 * @param res         The result uniformly drawn \ZnX polynomial.
 * @param nn          Number of coeffients in the polynomial (eq. degree of the cyclotomial poly)
 * @param low_bound   Lower bound of the generated numbers
 * @param high_bound  Upperbound of the generated numbers
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int uniform_random_pol_znx(PolyUniv* res, uint64_t nn, int64_t low_bound, int64_t high_bound);

/**
 * @brief Generates a random vector following a uniform distribution in res.
 *
 * In other words, every res_sl words/elements (of 64 bits), it fills limb_len
 * elements with a uniformly sampled random integer of nb_bits, between
 * [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param limb_len The size of each limb. A common value could be \N or \f$kN\f$ (to fill the A's)
 * @param res      The result.
 * @param nb_limbs The number of limbs. In other words, how many times we fill limb_len words.
 * @param res_sl   The stride (in elements) between two consecutive result limbs.
 *                 It indicates how many elements must be skipped in memory to reach
 *                 the start of the next limb in `res`.
 * @param nb_bits  The exponent of the distribution range = [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @retval -1 if an error occurs.
 * @retval  0 otherwise.
 *
 */
int uniform_random_vec(uint64_t limb_len, int64_t* res, uint64_t nb_limbs, uint64_t res_sl, uint64_t nb_bits);

/**
 * @brief Generates a random vector following a uniform distribution and return it in the DFT domain.
 *
 * @param module   Additionnal information for backend.
 * @param result_dft  The result in the DFT domain.
 * @param vec_size The number of limbs.
 * @param nb_bits  The exponent of the range = [-2^nb_bits, 2^nb_bits).
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 *
 */
int uniform_random_vec_znx_dft(const PvdaBackend* module, VecUnivDFT* result_dft, uint64_t vec_size, uint64_t nb_bits);

/**
 * @brief Generates a Random Vector following a normal distribution.
 *
 * @param res       The result.
 * @param res_size  The number of elements in the vector.
 * @param mu        The mean value of the distribution.
 * @param sigma     The standard deviation of the distribution.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.

 */
int normal_random_vec(double* res, uint64_t res_size, double mu, double sigma);

/**
 * @brief Adds a random normally sampled number to each element of vec
 *
 * res and vec can be the same vector for in-place addition
 *
 * @param res       The result.
 * @param vec_size  The number of elements in the vectors.
 * @param vec       The input vector.
 * @param mu        The mean parameter of the distribution.
 * @param sigma     The standard deviation of the distribution.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 *
 */
int add_normal_random_vec(double* res, size_t vec_size, const double* vec, double mu, double sigma);

#endif  // RNG_H
