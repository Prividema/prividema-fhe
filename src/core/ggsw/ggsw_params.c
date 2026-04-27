#include "ggsw_params.h"

#include <assert.h>
#include <stdlib.h>

#include "utils.h"

GGSWParams* new_ggsw_params(const GLWEParams* params_glwe, uint64_t k_tilde, uint64_t kappa_tilde,
                            uint64_t nb_limbs_tilde)
{
	assert(params_glwe);
	GGSWParams* params_ggsw = malloc(sizeof(GGSWParams));
	CHECK_ALLOC(params_ggsw, "malloc in new_ggsw_ct_params");

	params_ggsw->params_glwe               = params_glwe;
	params_ggsw->k_tilde                   = k_tilde;
	params_ggsw->kappa_tilde               = kappa_tilde;
	params_ggsw->ciphertext_nb_limbs_tilde = nb_limbs_tilde;

	return params_ggsw;
cleanup:
	return NULL;
}

void delete_ggsw_params(GGSWParams* params) { free(params); }

GLWEGadgetParams* new_glwegadget_params(const GLWEParams* params, uint64_t kappa_tilde, uint64_t l_tilde)
{
	assert(params);
	GLWEGadgetParams* params_glwegadget = malloc(sizeof(GLWEGadgetParams));
	CHECK_ALLOC(params_glwegadget, "malloc failed in new_ggsw_ct_params");

	params_glwegadget->params_glwe = params;
	params_glwegadget->kappa_tilde = kappa_tilde;
	params_glwegadget->l_tilde     = l_tilde;

	return params_glwegadget;
cleanup:
	return NULL;
}

void* delete_glwegadget_params(GLWEGadgetParams* params) { free(params); }

uint64_t glwegadget_coef_number(const GLWEGadgetParams* params_glwegadget)
{
	return params_glwegadget->l_tilde * glwe_coef_number(params_glwegadget->params_glwe);
}

uint64_t ggsw_num_rows(const GGSWParams* params) { return params->ciphertext_nb_limbs_tilde; };

uint64_t ggsw_params_l_tilde_a(const GGSWParams* params)
{
	return (params->ciphertext_nb_limbs_tilde + 1) / (params->k_tilde + 1);
}
uint64_t ggsw_params_l_tilde_b(const GGSWParams* params)
{
	return params->ciphertext_nb_limbs_tilde - params->k_tilde * ggsw_params_l_tilde_a(params);
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
