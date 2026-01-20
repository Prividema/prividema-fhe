#include "ggsw_ct_params.h"
#include <stdio.h>

/**
 * @brief Create a set of GGSW parameters.
 * 
 * @param params_glwe The GLWE parameters.
 * @param k_tilde The number of ZnX polynomial in the secret key.
 * @param kappa_tilde The 2-exponent of the base Bg_tilde.
 * @param n_limbs_tilde (k_tilde + 1)*l_tilde.
 * @return GGSWCtParams* 
 */
GGSWCtParams* new_ggsw_ct_params(GLWECtParams* params_glwe, uint64_t k_tilde, uint64_t kappa_tilde,  uint64_t n_limbs_tilde
){
    GGSWCtParams* params_ggsw = malloc(sizeof(GGSWCtParams));
    if (params_ggsw == NULL){
        perror("Malloc failed.");
        return NULL;
    }
    
    params_ggsw->params = params_glwe;
    params_ggsw->k_tilde = k_tilde;
    params_ggsw->kappa_tilde = kappa_tilde;
    params_ggsw->n_limbs_tilde = n_limbs_tilde;
    
    return params_ggsw;
}

/**
 * @brief Delete a GGSW parameters
 * 
 * @param params 
 */
void delete_ggsw_ct_params(GGSWCtParams* params){
    free(params);
}

/**
 * @brief Return the number of partialGGSW in a GGSW.
 * 
 * @param params The GGSW parameters.
 * @return int64_t 
 */
int64_t nb_partials(GGSWCtParams* params){
    return params->n_limbs_tilde/(params->k_tilde + 1);
}

/**
 * @brief Return the number of rows in a partialGGSW.
 * 
 * @param params The GGSW parameters.
 * @return int64_t 
 */
int64_t nb_rows_per_partial(GGSWCtParams* params){
    return params->k_tilde + 1;
}