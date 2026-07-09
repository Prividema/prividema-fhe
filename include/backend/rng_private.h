#ifndef RNG_PRIVATE_H
#define RNG_PRIVATE_H

#include "backend.h"

// =============================================
// |                                           |
// |              Random Numbers               |
// |                                           |
// =============================================

/**
 * @brief Generates a uniformly sampled random number in [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param module  The backend object
 * @param result  The resulting uniformly sampled integer
 * @param nb_bits The number of bits of the result
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ref_rand_uniform_pow2(const PvdaBackend* module, int64_t* result, uint64_t nb_bits);

/**
 * @brief Generates a uniformly sampled random number in [limit_down, limit_up] via
 * power-of-2 sampling and resampling if out-of-bounds
 *
 * @param module      The backend object
 * @param result      The resulting uniformly sampled integer
 * @param limit_down  The lower bound of the uniform sample
 * @param limit_up    The upper bound of the uniform sample
 *
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ref_rand_uniform(const PvdaBackend* module, int64_t* result, int64_t limit_down, int64_t limit_up);

/**
 * @brief Generates a uniformly random vector of numbers
 *
 * Coefficients are uniformly sampled in range [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param module  The backend object
 * @param res     The result uniformly drawn vector of numbers
 * @param nn      Number of coeffients in the polynomial (eq. degree of the cyclotomial poly)
 * @param nb_bits Number of randomness bits per coefficient.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int ref_rand_uniform_pow2_vec(const PvdaBackend* module, int64_t* res, uint64_t n, uint64_t nb_bits);

/**
 * @brief Generates a uniformly random binary vector of numbers
 *
 * Coefficients are uniformly sampled in {0, 1}
 *
 * @param module  The backend object
 * @param res     The result uniformly drawn \ZnX polynomial.
 * @param nn      Number of coeffients in the polynomial (eq. degree of the cyclotomial poly)
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int ref_rand_uniform_binary_vec(const PvdaBackend* module, uint64_t* res, uint64_t n);
#endif
