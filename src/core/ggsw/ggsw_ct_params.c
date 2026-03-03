#include "ggsw_ct_params.h"

#include <stdio.h>

#include "logger.h"

GGSWCtParams* new_ggsw_ct_params(const GLWECtParams* params_glwe, uint64_t k_tilde, uint64_t kappa_tilde,
                                 uint64_t n_limbs_tilde)
{
	GGSWCtParams* params_ggsw = malloc(sizeof(GGSWCtParams));
	if (log_is_null(params_ggsw, "malloc in new_ggsw_ct_params") < 0) return NULL;

	params_ggsw->params_glwe   = params_glwe;
	params_ggsw->k_tilde       = k_tilde;
	params_ggsw->kappa_tilde   = kappa_tilde;
	params_ggsw->n_limbs_tilde = n_limbs_tilde;

	return params_ggsw;
}

void delete_ggsw_ct_params(GGSWCtParams* params) { free(params); }

uint64_t nb_partials(const GGSWCtParams* params) { return params->n_limbs_tilde / (params->k_tilde + 1); }

uint64_t nb_rows_per_partial(const GGSWCtParams* params) { return params->k_tilde + 1; }