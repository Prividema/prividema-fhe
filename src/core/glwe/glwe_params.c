#include "glwe_params.h"

#include <stdlib.h>

#include "logger.h"
#include "utils.h"

GLWEParams* new_glwe_ct_params(uint64_t N, uint64_t k, uint64_t kappa, uint64_t n_limbs, double sigma)
{
	GLWEParams* params = malloc(sizeof(GLWEParams));
	CHECK_ALLOC(params, "params' malloc failed in new_glwe_ct_params");

	params->N       = N;
	params->k       = k;
	params->kappa   = kappa;
	params->n_limbs = n_limbs;
	params->sigma   = sigma;

	return params;
cleanup:
	return NULL;
}

void delete_glwe_ct_params(GLWEParams* params) { free(params); }
