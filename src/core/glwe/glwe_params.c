#include "glwe_params.h"

#include <stdlib.h>

#include "logger.h"
#include "utils.h"

GLWEParams* new_glwe_params(uint64_t nn, uint64_t k, uint64_t kappa, uint64_t n_limbs, double sigma)
{
	GLWEParams* params = malloc(sizeof(GLWEParams));
	CHECK_ALLOC(params, "params' malloc failed in new_glwe_ct_params");

	params->nn      = nn;
	params->k       = k;
	params->kappa   = kappa;
	params->n_limbs = n_limbs;
	params->sigma   = sigma;

	return params;
cleanup:
	return NULL;
}

void delete_glwe_params(GLWEParams* params) { free(params); }

uint64_t glwe_params_l(const GLWEParams* params_glwe) { return params_glwe->n_limbs / (params_glwe->k + 1); }

uint64_t glwe_params_n_limbs(const GLWEParams* params_glwe) { return params_glwe->n_limbs; }
