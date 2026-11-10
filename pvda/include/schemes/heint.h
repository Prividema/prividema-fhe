#ifndef HEINT_H
#define HEINT_H

#include <stdint.h>

#include "backend.h"
#include "glwe_params.h"
#include "univariate_polynomial.h"
#include "bivariate_polynomial.h"

/**
 * @brief Encodes "n_slots" integer slots in Z_t into the coefficient bivariate representation needed
 *        for HEInt arithmetic
 *
 * @param module The computation backend
 * @param params_glwe The GLWEParams, used for knowledge of N the polynomial number of coefficients and the depth of
 *               the bivariate output
 * @param out    The output bivariate polynomial
 * @param n_slots The number of elements to encode. For now only n_slots = N is supported
 * @param t      The size/modulus of the input ring Z_t
 * @param in     A vector of n_slots integers in Z_t
 *
 * @note For now only n_slots = N is supported
 * @note t must be a prime with t = 1 mod 2N and t <= 2^31, and ntt_table_generate(backend, t) must have been
 *       called beforehand, otherwise the function fails with "root table not found"
 *
 * @retval 0 if it succeeds
 * @retval -1 in case of error/failure
 *
 */
int heint_encode(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* out, uint64_t n_slots, int64_t t,
                 uint64_t* in);

/**
 * @brief Decodes integers from a bivariate polynomial into the "n_slots" slots in Z_t it represents
 *
 * @param module The computation backend
 * @param params_glwe The GLWEParams, used for knowledge of N the polynomial number of coefficients and the depth of
 *               the bivariate output
 * @param out    The output n_slots-length vector of elements of Z_t
 * @param n_slots The number of slots that were encoded. For now only n_slots = N is supported
 * @param t      The size/modulus of the input ring Z_t
 * @param in     A bivariate polynomial encoding some HEInt numbers
 *
 * @note For now only n_slots = N is supported
 * @note t must be a prime with t = 1 mod 2N and t <= 2^31, and ntt_table_generate(backend, t) must have been
 *       called beforehand, otherwise the function fails with "root table not found"
 *
 * @retval 0 if it succeeds
 * @retval -1 in case of error/failure
 *
 */
int heint_decode(const PvdaBackend* module, const GLWEParams* params_glwe, uint64_t* out, uint64_t n_slots, int64_t t,
                 PolyBiv* in);

int internal_slow_intt_heint(const PvdaBackend* module, uint64_t* out_int, uint64_t* in, uint64_t t);
int internal_slow_ntt_heint(const PvdaBackend* module, uint64_t* out_int, uint64_t* in, uint64_t t);
#endif
