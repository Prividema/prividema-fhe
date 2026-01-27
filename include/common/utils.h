#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>

/**
 * Generates an uniform random number.
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

#endif  // UTILS_H