#include <stdlib.h>
#include "glwe_ct_params.h"

GLWECtParams* new_glwe_ct_params(uint64_t N, uint64_t k, uint64_t kappa, uint64_t n_limbs
){
    GLWECtParams* params = malloc(sizeof(GLWECtParams));
    params->N = N;
    params->k = k;
    params->kappa = kappa;
    params->n_limbs = n_limbs;

    return params;
}

void delete_glwe_ct_params(GLWECtParams* params){
    free(params);
}
