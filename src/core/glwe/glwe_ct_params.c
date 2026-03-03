#include "glwe_ct_params.h"

#include <stdlib.h>

GLWECtParams* new_glwe_ct_params(uint64_t N, uint64_t k, uint64_t kappa, uint64_t n_limbs, double sigma)
{
	GLWECtParams* params = malloc(sizeof(GLWECtParams));
	params->N            = N;
	params->k            = k;
	params->kappa        = kappa;
	params->n_limbs      = n_limbs;
	params->sigma        = sigma;

	return params;
}

void delete_glwe_ct_params(GLWECtParams* params) { free(params); }
