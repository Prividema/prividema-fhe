#ifndef bivGGSW_CT_PARAMS_H
#define bivGGSW_CT_PARAMS_H

/**
 * @file ggsw_ct_params.h
 * @brief bivGGSW parameters
 *
 * In this header, we define the structure representing bivGGSW parameters.
 */

#include "glwe_ct_params.h"

/**
 * @struct GGSWCtParams
 * @brief bivGGSW parameters.
 */
typedef struct ggsw_ct_params
{
	const GLWECtParams* params_glwe;
	uint64_t k_tilde;  //!< k_tilde = 1 for RGSW (by default k = k_tilde=1). */
	uint64_t kappa_tilde;
	uint64_t n_limbs_tilde;  //!< n_limbs_tilde = (k_tilde + 1)*l_tilde. */
} GGSWCtParams;

/**
 * @brief Create a set of bivGGSW parameters.
 *
 * @param params_glwe The bivGLWE parameters.
 * @param k_tilde The number of Zn[X] polynomial in the secret key.
 * @param kappa_tilde The 2-exponent of the base (2^kappa_tilde).
 * @param n_limbs_tilde (k_tilde + 1)*l_tilde.
 *
 * @return GGSWCtParams*
 */
GGSWCtParams* new_ggsw_ct_params(const GLWECtParams* params, uint64_t k_tilde, uint64_t kappa_tilde, uint64_t n_limbs_tilde);

/**
 * @brief Delete a bivGGSW parameters
 *
 * @param params
 */
void delete_ggsw_ct_params(GGSWCtParams* params);

/**
 * @brief Return the number of partialGGSW in a bivGGSW.
 *
 * @param params The bivGGSW parameters.
 * @return int64_t
 */
uint64_t nb_partials(const GGSWCtParams* params);

/**
 * @brief Return the number of rows in a partialGGSW.
 *
 * @param params The bivGGSW parameters.
 * @return int64_t
 */
uint64_t nb_rows_per_partial(const GGSWCtParams* params);

/**
 * @struct PartialGGSWCtParams
 * @brief
 *
 */
typedef struct partialggsw_ct_params
{
	GLWECtParams* params_glwe;  //!< bivGLWE parameters. */
	uint64_t kappa_tilde;       //!< (2^kappa) = 2^-kappa_tilde. */
	uint64_t l_tilde;           //!< The precision. */
} PartialGGSWCtParams;

#endif  // bivGGSW_CT_PARAMS_H
