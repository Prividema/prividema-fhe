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
	uint64_t n_limbs_tilde;  //!< n_limbs_tilde = (k_tilde + 1)*l_tilde */

} GGSWCtParams;

int new_ggsw_ct_params(GGSWCtParams* res, GLWECtParams* params, uint64_t k_tilde, uint64_t kappa_tilde, uint64_t n_limbs_tilde);
void delete_ggsw_ct_params(GGSWCtParams* ct);
int64_t ggsw_size(GGSWCtParams* params);
int64_t nb_partials(GGSWCtParams* params);
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
