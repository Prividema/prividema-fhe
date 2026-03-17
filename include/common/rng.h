#ifndef RNG_H
#define RNG_H

#include <stdint.h>

#include "maths_structures.h"
#include "spqlios_alias.h"

// =============================================
// |                                           |
// |              Random Numbers               |
// |                                           |
// =============================================

/**
 * @brief Generates an uniform random number.
 *
 * @param result  A pointer that Point to the generated value.
 * @param nb_bits The number of bits of the result, sign bit included.
 *
 * @retval • `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval • `0` otherwise.
 */
int rand_uniform(int64_t* result, uint64_t nb_bits);

/**
 * @brief Generates a gaussian random number with the given parameters.
 *
 * @param result A pointer that Point to the generated value.
 * @param mu     The mean value.
 * @param sigma  The standard deviation.
 *
 * @retval • `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval • `0` otherwise.
 */
int rand_normal(double* result, double mu, double sigma);

// =============================================
// |                                           |
// |      Random Polynomials and Vectors       |
// |                                           |
// =============================================

/**
 * @brief Generates a uniformly random Zn[X] polynomial.
 *
 * @param res     The result uniformly drawn Zn[X] polynomial.
 * @param N       The degree of the cyclotomic polynomial.
 * @param nb_bits The exponent of the distribution range = [-2^nb_bits, 2^nb_bits).
 *
 * @retval • `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval • `0` otherwise.
 */
int uniform_random_pol_znx(PolyUniv* res, uint64_t N, uint64_t nb_bits);

/**
 * @brief Generates a random vector following a uniform distribution in res.
 *
 * @param limb_len The size of each limb. N for TGLWE and TGGSW.
 * @param res      The result.
 * @param nb_limbs The number of limbs.
 * @param res_sl   The stride (in elements) between two consecutive result limbs.
 *                 It indicates how many elements must be skipped in memory to reach
 *                 the start of the next limb in `res`.
 * @param nb_bits  The exponent of the distribution range = [-2^nb_bits, 2^nb_bits).
 *
 * @retval • `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval • `0` otherwise.
 *
 * @note For each call to this function, we'll generate a new seed.
 */
int uniform_random_vec(uint64_t limb_len, int64_t* res, int64_t nb_limbs, int64_t res_sl, uint64_t nb_bits);

/**
 * @brief Generates a random vector following a uniform distribution and return it in the DFT domain.
 *
 * @param module   Additionnal information for backend.
 * @param result_dft  The result in the DFT domain.
 * @param vec_size The number of limbs.
 * @param nb_bits  The exponent of the range = [-2^nb_bits, 2^nb_bits).
 *
 * @retval • `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval • `0` otherwise.
 *
 * @note For each call to this function, we'll generate a new seed.
 */
int uniform_random_vec_znx_dft(const MODULE* module, VecUnivDFT* result_dft, uint64_t vec_size, uint64_t nb_bits);

/**
 * @brief Generates a Random Vector following a normal distribution.
 *
 * @param res       The result.
 * @param res_size  The number of elements in the vector.
 * @param mu        The mean value of the distribution.
 * @param sigma     The standard deviation of the distribution.
 *
 * @retval • `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval • `0` otherwise.
 *
 * @note For each call to this function, we'll generate a new seed.
 */
int normal_random_vec(double* res, int64_t res_size, double mu, double sigma);

#endif  // RNG_H
