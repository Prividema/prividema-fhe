#include <stdlib.h>
#include "glwe_ct_params.h"

/**
 * @brief Creates a new set of GLWE parameters.
 * 
 * @param N The polynomial maximum degree in X
 * @param k The number of ZnX polynomial in the secret key.
 * @param kappa The exponent for the base-2^kappa representation.
 * @param n_limbs (k+1)*l.
 * @param sigma The standard deviation of the error distribution.
 * @return GLWECtParams* 
 */
GLWECtParams* new_glwe_ct_params(uint64_t N, uint64_t k, uint64_t kappa, uint64_t n_limbs, uint64_t sigma
){
    GLWECtParams* params = malloc(sizeof(GLWECtParams));
    params->N = N;
    params->k = k;
    params->kappa = kappa;
    params->n_limbs = n_limbs;
    params->sigma = sigma;

    return params;
}

/**
 * @brief Deletes the set of GLWE parameters.
 * 
 * @param params The GLWE parameters.
 */
void delete_glwe_ct_params(GLWECtParams* params){
    free(params);
}
