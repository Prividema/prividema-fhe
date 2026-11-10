#ifndef MONTARITH_H
#define MONTARITH_H

#include <stdint.h>

/**
 * @brief Encodes a number x to its montgomery form for radix of 32bits (R = 2^32)
 *
 * q being the modulus used, it must be at most 2^31
 *
 * @param x       The number to encode
 * @param q       The modulus to encode for, N in the original 1985 Montgomery paper
 * @param q_tild  -q^-1 mod R, that is, N' in the original paper
 *                One would obtain said value using montgomery_tild_32bit
 * @param r2modq  (R^2 mod q), obtained using montgomery_r2modq_32bit
 *
 * @info    This function DOES NOT use and should not use xR^-1 mod q.
 *          Said operation is both slow and non-constant time depending on input x.
 *          Input x at some points will probably be cleartext/plaintext data.
 *          Thus, this implementation uses r2modq for a CONSTANT-TIME encoding
 *          to avoid a side-channel.
 *
 * @returns x in montgomery form
 */
uint64_t montgomery_encode_32bit(uint64_t x, uint64_t q, uint64_t q_tild, uint64_t r2modq);

/**
 * @brief Computes -q^-1 mod R (with R = 2^32) (N' in the original 1985 Montgomery paper)
 *
 * @param q modulus
 *
 * @return q_tild such that q_tild*q mod R is congruent to -1
 */
uint64_t montgomery_tild_32bit(uint64_t q);

/**
 * @brief Computes (R^2 mod q), which is used for constant-time montgomery encoding, for R = 2^32
 *
 * This function is not constant-time on machines that have non-constant integer modulus
 * operation.
 * q should in most cases be a public parameter, so this should not be an issue.
 *
 * @param q The used modulus
 *
 * @returns (R^2 mod q)
 */
uint64_t montgomery_r2modq_32bit(uint64_t q);

/**
 * @brief Implements operation REDC from the montgomery paper, with radix R = 2^32
 *
 * @param m       The input value
 * @param q       The modulus, N in the paper
 * @param q_tild  N' in the paper, obtained by montgomery_tild_32bit
 *
 * @return REDC(m) for modulus q
 */
uint64_t montgomery_red_32bit(uint64_t m, uint64_t q, uint64_t q_tild);

/**
 * @brief Decodes a number in montgomery form back to its normal representation
 *
 * @param x A number in montgomery form
 * @param q       The modulus, N in the paper
 * @param q_tild  N' in the paper, obtained by montgomery_tild_32bit
 *
 * @returns The number that x represents, in normal form (not Montgomery)
 */
uint64_t montgomery_decode_32bit(uint64_t x, uint64_t q, uint64_t q_tild);

/**
 * @brief Multiplies two numbers in montgomery form and returns the result in montgomery form
 *
 * @param x_m A number in montgomery form, representing x
 * @param y_m A number in montgomery form, representing y
 * @param q       The modulus, N in the paper
 * @param q_tild  N' in the paper, obtained by montgomery_tild_32bit
 *
 * @returns The product of x and y in montgomery form
 */
uint64_t montgomery_mult_32bit(uint64_t x_m, uint64_t y_m, uint64_t q, uint64_t q_tild);

/**
 * @brief Computes the power of a number in montgomery form
 *
 * @param base_m  the base, in montgomery form
 * @param exp     the exponent, in normal form, that is, NOT Montgomery
 * @param q       The modulus, N in the paper
 * @param q_tild  N' in the paper, obtained by montgomery_tild_32bit
 * @param one_m   The number 1 in montgomery form for the provided modulus
 *
 * @returns base^exp in montgomery form
 *
 * @warning This function uses the square-and-multiply algorithm and thus is not constant time.
 */
uint64_t montgomery_pow_exp_32bit(uint64_t base_m, uint64_t exp, uint64_t q, uint64_t q_tild, uint64_t one_m);

#endif
