#ifndef RNG_H
#define RNG_H

#include <stdint.h>
#include "maths_structures.h"
#include "spqlios_alias.h"

// Random Numbers

/**
 * Generates a uniform random number.
 *
 * @param result A pointer that will point to the generated value.
 * @param nb_bits The number of bits of the result, sign bit included.
 * 
 * @retval `-1` if an error occurs. A message will be printed in this case.
 * @retval `0` otherwise.
 */
int rand_uniform(int64_t *result, uint64_t nb_bits);

/**
 * Generates a gaussian random number with the given parameters.
 *
 * @param result A pointer that will point to the generated value.
 * @param mu     The mean value.
 * @param sigma  The standard deviation.
 *
 * @retval `-1` if an error occurs. A message will be printed in this case.
 * @retval `0` otherwise.
 */
int rand_normal(double *result, double mu, double sigma);


// Random Vectors

/**
 * @brief Generates a random vector following a uniform distribution in res.
 * 
 * @param limb_len The size of each limb. N for TGLWE and TGGSW.
 * @param res      The result.
 * @param res_size The number of limbs.
 * @param res_sl   The stride (in elements) between two consecutive result limbs.
 *                 It indicates how many elements must be skipped in memory to reach
 *                 the start of the next limb in `res`.
 * @param nb_bits  The exponent of the distribution range = [-2^nb_bits, 2^nb_bits).
 *  
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int inplace_uniform_random_vec(uint64_t limb_len, int64_t* res, int64_t limb_nb, int64_t res_sl, uint64_t nb_bits);

/**
 * @brief Generates a random vector - with coef_nb int coefficients - following a uniform distribution U([-2^nb_bits, 2^nb_bits)).
 * 
 * @param coef_nb The number of coefficient of the result vector.
 * @param nb_bits The exponent of the distribution range = [-2^nb_bits, 2^nb_bits).
 *  
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int64_t* new_uniform_random_vec(uint64_t coef_nb, uint64_t nb_bits);

/**
 * @brief Generates a random vector following a uniform distribution and return it in DFT space.
 * 
 * @param module   The module holding the degree N in X.
 * @param res_dft  The result in DFT space.
 * @param res_size The number of limbs.
 * @param nb_bits  The exponent of the range = [-2^nb_bits, 2^nb_bits).
 *  
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int inplace_uniform_random_vec_znx_dft(MODULE* module, PolyUnivDFT* res_dft, int64_t res_size, uint64_t nb_bits);

/**
 * @brief Generates a random Zn[X] vector, following a uniform distribution and returns it in DFT space .
 * 
 * @param module 
 * @param vec_size 
 * @param nb_bits 
 * @return VecUnivDFT* 
 */
VecUnivDFT* new_uniform_random_vec_znx_dft(MODULE* module, uint64_t vec_size, uint64_t nb_bits);

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

#endif // RNG_H