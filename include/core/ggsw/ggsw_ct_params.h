#ifndef GGSW_CT_PARAMS_H
#define GGSW_CT_PARAMS_H

/**
 * @file ggsw_ct_params.h
 * @brief GGSW parameters
 * 
 * In this header, we define the structure representing GGSW parameters.
 */

#include "glwe_ct_params.h"

/**
 * @struct GGSWCtParams 
 * @brief GGSW parameters.
 */
typedef struct ggsw_ct_params {
	GLWECtParams* params;    //!< GLWE parameters. */
	uint64_t k_tilde;        //!< k_tilde = 1 for RGSW (by default k=k_tilde=1). */ 
	uint64_t kappa_tilde;    //!< Bg = 2^-kappa_tilde. */
	uint64_t n_limbs_tilde;  //!< n_limbs_tilde = (k_tilde + 1)*l_tilde. */
} GGSWCtParams;

/**
 * @brief Create a set of GGSW parameters.
 * 
 * @param params_glwe The GLWE parameters.
 * @param k_tilde The number of ZnX polynomial in the secret key.
 * @param kappa_tilde The 2-exponent of the base Bg_tilde.
 * @param n_limbs_tilde (k_tilde + 1)*l_tilde.
 * 
 * @return GGSWCtParams* 
 */
GGSWCtParams* new_ggsw_ct_params(GLWECtParams* params,
								 uint64_t k_tilde,
								 uint64_t kappa_tilde,
								 uint64_t n_limbs_tilde);

/**
 * @brief Delete a GGSW parameters
 * 
 * @param params 
 */
void delete_ggsw_ct_params(GGSWCtParams* params);

/**
 * @brief Return the number of partialGGSW in a GGSW.
 * 
 * @param params The GGSW parameters.
 * @return int64_t 
 */
int64_t nb_partials(GGSWCtParams* params);

/**
 * @brief Return the number of rows in a partialGGSW.
 * 
 * @param params The GGSW parameters.
 * @return int64_t 
 */
int64_t nb_rows_per_partial(GGSWCtParams* params);

/**
 * @struct PartialGGSWCtParams
 * @brief 
 * 
 */
typedef struct partialggsw_ct_params {
	GLWECtParams* params;	//!< GLWE parameters. */
	uint64_t kappa_tilde;   //!< Bg = 2^-kappa_tilde. */
	uint64_t l_tilde; 		//!< The precision. */
} PartialGGSWCtParams;

#endif  // GGSW_CT_PARAMS_H
