#ifndef MATH_UTILS
#define MATH_UTILS

#include <stdint.h>

/**
 * @brief Computes the inverse of x in Z_mod.
 * Undefined behaviour if x is not coprime with mod
 *
 * @param x an input value coprime with mod
 * @param mod The modulus of the ring
 *
 * @returns y such that xy is congruent to 1 modulo mod
 *
 */
uint64_t mod_inv(int64_t x, int64_t mod);

/**
 * @brief Computes the distance between two torus elements
 *
 * In other words, abs( (a mod 1) - (b mod 1))
 * Mostly used in tests
 *
 * @param a An element in the torus
 * @param b Another element in the torus
 *
 * @returns The (shortest) distance between the two elements in the torus
 */
double rnx_torus_distance(double a, double b);

/**
 * @brief Computes the distance between two torus elements
 *
 * In other words, abs( (a mod 1) - (b mod 1))
 * Mostly used in tests
 *
 * @param a An element in the torus
 * @param b Another element in the torus
 *
 * @returns The (shortest) distance between the two elements in the torus
 */
uint64_t tnx_torus_distance(uint64_t a, uint64_t b);

/**
 * @brief uint64_t implementation for INT_ROUND_UP_DIV
 *
 */
uint64_t u64_round_up_div(uint64_t num, uint64_t step);

/**
 * @brief int32_t implementation for INT_ROUND_UP_DIV
 *
 */
int32_t i32_round_up_div(int32_t num, int32_t step);

/**
 * @brief Generic macro to compute an integer division with rounding up when the remainder is non-zero.
 *
 * @param num The dividend
 * @param step The divisior
 *
 */
#define INT_ROUND_UP_DIV(num, step) _Generic((num), uint64_t: u64_round_up_div, int32_t: i32_round_up_div)(num, step)

/**
 * @brief Returns ceil(log2(num)), with a 0 return value if num == 0.
 *
 * @param num The input number
 * @returns ceil(log_2(num)) or 0 if num == 0
 */
uint64_t next_pow2_log(uint64_t num);

#endif
