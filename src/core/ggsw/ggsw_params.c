#include "ggsw_params.h"

#include <assert.h>
#include <stdio.h>

#include "logger.h"
#include "utils.h"

GGSWParams* new_ggsw_params(const GLWEParams* params_glwe, uint64_t k_tilde, uint64_t kappa_tilde,
                            uint64_t n_limbs_tilde)
{
	assert(params_glwe);
	GGSWParams* params_ggsw = malloc(sizeof(GGSWParams));
	CHECK_ALLOC(params_ggsw, "malloc in new_ggsw_ct_params");

	params_ggsw->params_glwe   = params_glwe;
	params_ggsw->k_tilde       = k_tilde;
	params_ggsw->kappa_tilde   = kappa_tilde;
	params_ggsw->n_limbs_tilde = n_limbs_tilde;

	return params_ggsw;
cleanup:
	return NULL;
}

void delete_ggsw_params(GGSWParams* params) { free(params); }

uint64_t ggsw_num_pggsw(const GGSWParams* params) { return params->n_limbs_tilde / (params->k_tilde + 1); }

uint64_t ggsw_num_rows_per_pggsw(const GGSWParams* params) { return params->k_tilde + 1; }
