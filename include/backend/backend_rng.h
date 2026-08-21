#ifndef PVDA_BACK_RNG_H
#define PVDA_BACK_RNG_H

#include "backend.h"

/**
 * @brief Generates a uniformly sampled random number in [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param result  The resulting uniformly sampled integer
 * @param nb_bits The number of bits of the result
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int pvda_rand_uniform_pow2(const PvdaBackend* backend, int64_t* result, uint64_t nb_bits);

/**
 * @brief Generates a uniformly sampled vector of random numbers in [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param result  The resulting uniformly sampled integer
 * @param n       Number of elements in the vector
 * @param nb_bits The number of bits of the result
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int pvda_rand_uniform_pow2_vec(const PvdaBackend* backend, int64_t* result, uint64_t n, uint64_t nb_bits);

/**
 * @brief Generates a uniformly sampled vector of random numbers in {0,1}
 *
 * @param result  The resulting uniformly sampled integer
 * @param n       Number of elements in the vector
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int pvda_rand_uniform_binary_vec(const PvdaBackend* backend, uint64_t* result, uint64_t n);

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
int pvda_rand_uniform(const PvdaBackend* module, int64_t* result, int64_t limit_down, int64_t limit_up);

#endif
