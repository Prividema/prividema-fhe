#include <stdlib.h>
#include "glwe_ct_params.h"

int new_glwe_ct_params(GLWECtParams* res, uint64_t N, uint64_t k, uint64_t kappa, uint64_t n_limbs
){
    res = malloc(sizeof(GLWECtParams));
    if(res == NULL)
        return -1;
    res->N = N;
    res->k = k;
    res->kappa = kappa;
    res->n_limbs = n_limbs;

    return 0;
}

void delete_glwe_ct_params(GLWECtParams* params){
    free(params);
}
