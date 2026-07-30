#ifndef HEINT_H
#define HEINT_H

#include <stdint.h>

#include "backend.h"
#include "glwe_params.h"
#include "maths_structures.h"

/**
 * @brief Encodes "n_slots" integer slots in Z_t into the coefficient bivariate representation needed
 *        for HEInt arithmetic
 *
 * @param backend The computation backend
 * @param params The GLWEParams, used for knowlege of N the polynomial number of coefficients and the depth of
 *               the bivariate output
 * @param out    The output bivariate polynomial
 * @param n_slots The number of elements to encode. For now only n_slots = N is supported
 * @param t      The size/moduli of the input ring Z_t
 * @param in     A vector of n_slots integers in Z_t
 *
 * @note For now only n_slots = N is supported
 *
 * @retval 0 if it succeeds
 * @retval -1 in case of error/failure
 *
 */
int heint_encode(const PvdaBackend* backend, const GLWEParams* params, PolyBiv* out, uint64_t n_slots, int64_t t,
                 uint64_t* in);

/**
 * @brief Decodes integers from a bivariate polynomial into the "n_slots" slots in Z_t it represents
 *
 * @param backend The computation backend
 * @param params The GLWEParams, used for knowlege of N the polynomial number of coefficients and the depth of
 *               the bivariate output
 * @param out    The output n_slots-lenght vector of elements of Z_t
 * @param n_slots The number of that were encoded. For now only n_slots = N is supported
 * @param t      The size/moduli of the input ring Z_t
 * @param in     A bivariate polynomial encoding some HEInt numbers
 *
 * @note For now only n_slots = N is supported
 *
 * @retval 0 if it succeeds
 * @retval -1 in case of error/failure
 *
 */
int heint_decode(const PvdaBackend* backend, const GLWEParams* params, uint64_t* out, uint64_t n_slots, int64_t t,
                 PolyBiv* in);

int internal_slow_intt_heint(const PvdaBackend* backend, uint64_t* out_int, uint64_t* in, uint64_t t);
int internal_slow_ntt_heint(const PvdaBackend* backend, uint64_t* out_int, uint64_t* in, uint64_t t);
#endif
