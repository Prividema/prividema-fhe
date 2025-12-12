#include "distributions.h"
#include "utils.h"

/**
 * @brief Generates a Random Vector following an Uniform distribution.
 * 
 * @param limb_size The size of each limb. N for TGLWE and TGGSW.
 * @param res The result.
 * @param n_limbs The number of limbs.
 * @param res_sl The stride (in elements) between two consecutive result limbs.
 *               It indicates how many elements must be skipped in memory to reach
 *               the start of the next limb in `res`.
 *  
 * @retval - `-1` if an error occurs. In this case the error is from a syscall in rand_uniform.
 * @retval - `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int uniform_random_vec(int64_t limb_len, int64_t* res, int64_t n_limbs, int64_t res_sl)
{
    for(int i = 0; i < n_limbs; i++)
        for(int j = 0; j < limb_len ; j++){
            
            if(rand_uniform(res + i*res_sl + j) < 0) // TODO : def a range or mod by 2^kappa the result
                return -1;

            }

    return 0;
}

/**
 * @brief Generates a Random Vector following an Uniform distribution.
 * 
 * @param limb_size The size of each limb. N for TGLWE and TGGSW.
 * @param res The result.
 * @param n_limbs The number of limbs.
 * @param res_sl The stride (in elements) between two consecutive result limbs.
 *             It indicates how many elements must be skipped in memory to reach
 *             the start of the next limb in `res`.
 *  
 * @retval - `-1` if an error occurs. In this case the error is from a syscall in rand_uniform.
 * @retval - `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int normal_random_vec(int64_t limb_len, int64_t* res, int64_t n_limbs, int64_t res_sl)
{
    for(int i = 0; i < n_limbs; i++)
        for(int j = 0; j < limb_len ; j++){
            
            if(rand_normal((double*)(res + i*res_sl + j),0,1) < 0) // TODO : decompose to xy
                return -1;

            }

    return 0;
}