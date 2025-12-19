#include "ggsw_ct_params.h"

/**
 * @brief Create a GGSW parameters
 * 
 * @param params The 
 */
int new_ggsw_ct_params(GGSWCtParams* res,
                        GLWECtParams* params, 
                        uint64_t k_tilde, uint64_t kappa_tilde,  uint64_t n_limbs_tilde
){
    res = malloc(sizeof(GGSWCtParams));
    if(res == NULL)
        return -1;
    res->params = params;
    res->k_tilde = k_tilde;
    res->kappa_tilde = kappa_tilde;
    res->n_limbs_tilde = n_limbs_tilde;
    
    return 0;
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