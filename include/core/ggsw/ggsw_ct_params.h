#ifndef bivGGSW_CT_PARAMS_H
#define bivGGSW_CT_PARAMS_H

/**
 * @file ggsw_ct_params.h
 *
 * In this header, we define the structure representing bivGGSW parameters.
 */

#include "glwe_ct_params.h"

/**
 * @struct GGSWCtParams
 * 
 * @brief GGSW parameters.
 */
typedef struct ggsw_ct_params
{
	/// GLWE Parameters.
	const GLWECtParams* params_glwe;
	uint64_t k_tilde; 		// k_tilde = 1 for RGSW (by default k=k_tilde=1)
	uint64_t kappa_tilde;
	uint64_t n_limbs_tilde; // n_limbs_tilde = (k_tilde + 1) * l_tilde.
} GGSWCtParams;

/**
 * @brief Creates a set of GGSW Parameters.
 *
 * @param params_glwe 	The GLWE parameters.
 * @param k_tilde 	  	The number of polynomials in the secret key.
 * @param kappa_tilde 	The 2-exponent of the base (2^kappa_tilde).
 * @param n_limbs_tilde (k_tilde + 1) * l_tilde.
 *
 * @return A Pointer to the set of GGSW Parameters.
 */
GGSWCtParams* new_ggsw_ct_params(const GLWECtParams* params, uint64_t k_tilde, uint64_t kappa_tilde, uint64_t n_limbs_tilde);

/**
 * @brief Deletes a GGSW parameters.
 *
 * @param params A Pointer to the GGSW parameters.
 */
void delete_ggsw_ct_params(GGSWCtParams* params);

/**
 * @brief Gets the number of partialGGSW in a GGSW.
 *
 * @param params The GGSW parameters.
 * 
 * @return The number of partialGGSW in a GGSW.
 */
uint64_t nb_partials(const GGSWCtParams* params);

/**
 * @brief Gets the number of rows in a partialGGSW.
 *
 * @param params The GGSW parameters.
 * 
 * @return The number of rows in a partialGGSW.
 */
uint64_t nb_rows_per_partial(const GGSWCtParams* params);

/**
 * @struct PartialGGSWCtParams
 * 
 * @brief Partial GGSW Parameters.
 */
typedef struct partialggsw_ct_params
{
	/// GLWE parameters.
	GLWECtParams* params_glwe;  
	uint64_t kappa_tilde;       // (2^kappa) = 2^-kappa_tilde.
	uint64_t l_tilde;           // The precision.
} PartialGGSWCtParams;

#endif  // bivGGSW_CT_PARAMS_H
