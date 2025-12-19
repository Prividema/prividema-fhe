#include "ggsw_ciphertext.h"

/**
 * @brief Creates a bivGGSW, filled with 0.
 * 
 * @param res The result GGSW ciphertext.
 * @param params The GGSW parameters.
 */
int new_ggsw(GGSWCiphertext* res, GGSWCtParams* params, MatBiv* ct
){
    int64_t limb_size = params->params->N;
    int64_t nb_rows = params->n_limbs_tilde;
    int64_t nb_cols = params->n_limbs_tilde;

    res = malloc(sizeof(GGSWCiphertext));
    if(res == NULL)
        return -1;

    res->params = params;   

    if(ct == NULL)
        res->ct = calloc(nb_rows * nb_cols * limb_size, sizeof(int64_t));
        if(res->ct == NULL)
            return -1;
    else
        res->ct = ct;
    
    return 0;
}

/**
 * @brief Delete a bivGGSW ciphertext.
 * 
 * @param ct The bivGGSW ciphertext.
 */
void delete_ggsw(GGSWCiphertext* ct
){
    delete_ggsw_ct_params(ct->params);
    free(ct);
}


/**
 * @brief Creates a bivGGSW in DFT space, filled with 0.
 * 
 * @param res The result bivGGSW ciphertext in DFT space.
 * @param params The GGSW parameters.
 */
int new_ggsw_prepared(GGSWPreparedCt* res, GGSWCtParams* params, MatBivDFT* ct
){
    int64_t limb_size = params->params->N;
    int64_t nb_rows = params->n_limbs_tilde;
    int64_t nb_cols = params->n_limbs_tilde;

    res = malloc(sizeof(GGSWCiphertext));
    if(res == NULL)
        return -1;

    res->params = params;   

    if(ct == NULL)
        res->ct = calloc(nb_rows * nb_cols * limb_size, sizeof(double));
        if(res->ct == NULL)
            return -1;
    else
        res->ct = ct;
    
    return 0;
}

/**
 * @brief Delete a bivGGSW in DFT space.
 * 
 * @param ct The bivGGSW in DFT space.
 */
void delete_ggsw_prepared(GGSWPreparedCt* ct
){
    delete_ggsw_ct_params(ct->params);
    free(ct);
}