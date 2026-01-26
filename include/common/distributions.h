#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <stdint.h>
#include "utils.h"
#include "structure_alias.h"
#include "spqlios_alias.h"

/**
 * @brief Generates a random vector following a uniform distribution.
 * 
 * @param limb_len The size of each limb. N for TGLWE and TGGSW.
 * @param res      The result.
 * @param res_size The number of limbs.
 * @param res_sl   The stride (in elements) between two consecutive result limbs.
 *                 It indicates how many elements must be skipped in memory to reach
 *                 the start of the next limb in `res`.
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 *  
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int new_uniform_random_vec(uint64_t limb_len, int64_t* res, int64_t res_size, int64_t res_sl, uint64_t nb_bits);

/**
 * @brief Generates a random vector following a uniform distribution and return it in DFT space.
 * 
 * @param module   The module holding the degree N in X.
 * @param res_dft      The result in DFT space.
 * @param res_size The number of limbs.
 * @param nb_bits  The exponent of the range = [-2^nb_bits, 2^nb_bits).
 *  
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int new_uniform_random_vec_dft(MODULE* module, PolyUnivDFT* res_dft, int64_t res_size, uint64_t nb_bits);

/**
 * @brief Generates a Random Vector following a normal distribution.
 * 
 * @param limb_len The size of each limb. N for TGLWE and TGGSW.
 * @param res The result.
 * @param res_size The number of limbs.
 * @param res_sl The stride (in elements) between two consecutive result limbs.
 *             It indicates how many elements must be skipped in memory to reach
 *             the start of the next limb in `res`.
 * @param mu The mean value of the distribution.
 * @param sigma The variance of the distribution.
 *  
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int new_normal_random_vec(uint64_t limb_len, double* res, int64_t res_size, int64_t res_sl, double mu, double sigma);

#endif // DISTRIBUTION_H