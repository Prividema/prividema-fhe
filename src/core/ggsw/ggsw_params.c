#include "ggsw_params.h"

#include <assert.h>
#include <stdlib.h>

#include "utils.h"

GGSWParams* new_ggsw_params(const GLWEParams* params_glwe, uint64_t k_tilde, uint64_t nb_gglwes)
{
	assert(params_glwe);
	GGSWParams* params_ggsw = malloc(sizeof(GGSWParams));
	CHECK_ALLOC(params_ggsw, "malloc in new_ggsw_ct_params");

	params_ggsw->params_glwe          = params_glwe;
	params_ggsw->k_tilde              = k_tilde;
	params_ggsw->ciphertext_nb_gglwes = nb_gglwes;

	return params_ggsw;
cleanup:
	return NULL;
}

void delete_ggsw_params(GGSWParams* params) { free(params); }

uint64_t ggsw_num_rows(const GGSWParams* params) { return params->ciphertext_nb_gglwes; };

uint64_t ggsw_params_l_tilde_a(const GGSWParams* params)
{
	return (params->ciphertext_nb_gglwes + 1) / (params->k_tilde + 1);
}
uint64_t ggsw_params_l_tilde_b(const GGSWParams* params)
{
	return params->ciphertext_nb_gglwes - params->k_tilde * ggsw_params_l_tilde_a(params);
}

uint64_t ggsw_coef_number(const GGSWParams* params_ggsw)
{
	return ggsw_num_rows(params_ggsw) * glwe_coef_number(params_ggsw->params_glwe);
}

uint64_t ggsw_coef_number_dft(const GGSWParams* params_ggsw)
{
	return (ggsw_num_rows(params_ggsw) * glwe_coef_number(params_ggsw->params_glwe)) / 2;
}

uint64_t ggsw_bytes(const GGSWParams* params_ggsw)
{
	int64_t N = params_ggsw->params_glwe->nn;
	return ggsw_total_n_glwe_limbs(params_ggsw) * N * sizeof(int64_t);
}

uint64_t ggsw_total_n_glwe_limbs(const GGSWParams* params_ggsw)
{
	return ggsw_num_rows(params_ggsw) * glwe_params_n_limbs(params_ggsw->params_glwe);
}
