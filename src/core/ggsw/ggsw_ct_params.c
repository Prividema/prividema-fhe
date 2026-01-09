#include "ggsw_ct_params.h"
#include <errno.h>

/**
 * @brief Create a GGSW parameters
 * 
 * @param params The 
 */
GGSWCtParams* new_ggsw_ct_params(uint64_t k_tilde, uint64_t kappa_tilde,  uint64_t n_limbs_tilde
){
    GGSWCtParams* params = malloc(sizeof(GGSWCtParams));
    if (params == NULL)
        perror("Malloc failed.");
        return NULL;
    
    params->params = params;
    params->k_tilde = k_tilde;
    params->kappa_tilde = kappa_tilde;
    params->n_limbs_tilde = n_limbs_tilde;
    
    return params;
}

/**
 * @brief Delete a GGSW parameters
 * 
 * @param params 
 */
void delete_ggsw_ct_params(GGSWCtParams* params){
    delete_glwe_ct_params(params->params);
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