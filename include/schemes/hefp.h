#ifndef HEFP_H
#define HEFP_H

#include "univariate_polynomial.h"
#ifndef _cplusplus
#include <complex.h>
#endif
#include <stdint.h>

#include "backend.h"
#include "glwe_params.h"
#include "bivariate_polynomial.h"

/**
 * @brief Encodes "n_slots" complex numbers using scale 2^-scale_bits from its slots representation to a
 * coefficient representation.
 *
 * @param module The computation backend
 * @param params_glwe The GLWEParams
 * @param out    The output bivariate polynomial
 * @param n_slots The number of elements to encode. Should be a power of 2 and at most N/2.
 * @param scale_bits The scale factor \Delta is 2^-scale_bits
 * @param in A vector of n_slots complex values
 *
 * @retval 0 if it succeeds
 * @retval -1 if a memory allocation fails
 *
 */
int hefp_encode(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* out, uint64_t n_slots,
                int64_t scale_bits, double _Complex* in);

/**
 * @brief Decodes a bivariate polynomial encoding HEFixedPoint n_slots slots back into the n_slots complex values
 *
 * @param module The computation backend
 * @param params_glwe The GLWEParams
 * @param out    The output n_slots-length vector of complex values
 * @param n_slots The number of elements encoded. Should be a power of 2 and at most N/2.
 * @param scale_bits The scale factor \Delta is 2^-scale_bits
 * @param in A bivariate polynomial encoding some HEFixedPoint numbers
 *
 * @retval 0 if it succeeds
 * @retval -1 if a memory allocation fails
 *
 */
int hefp_decode(const PvdaBackend* module, const GLWEParams* params_glwe, double _Complex* out, uint64_t n_slots,
                int64_t scale_bits, PolyBiv* in);
#endif
