#include "distributions.h"
#include "utils.h"
/**
 * @brief Generates a Random Vector following an Uniform distribution.
 * 
 * @param limb_len The size of each limb. N for TGLWE and TGGSW.
 * @param res The result.
 * @param res_size The number of limbs.
 * @param res_sl The stride (in elements) between two consecutive result limbs.
 *               It indicates how many elements must be skipped in memory to reach
 *               the start of the next limb in `res`.
 * @param nb_bits The exponent of the range = [-2^e, 2^e].
 *  
 * @retval `-1` if an error occurs. In this case the error is from a syscall in rand_uniform.
 * @retval `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int uniform_random_vec(int64_t limb_len, 
                       int64_t* res, int64_t res_size, int64_t res_sl, 
                       int nb_bits)
{
    for(int i = 0; i < res_size; i++)
    {
        for(int j = 0; j < limb_len ; j++)
        {
            if(rand_uniform(res + i*res_sl + j, nb_bits) < 0) // TODO : def a range 
                return -1;
        }
    }
    return 0;
}

/**
 * @brief Generates a Random Vector following an Uniform distribution.
 * 
 * @param limb_len The size of each limb. N for TGLWE and TGGSW.
 * @param res The result.
 * @param res_size The number of limbs.
 * @param res_sl The stride (in elements) between two consecutive result limbs.
 *             It indicates how many elements must be skipped in memory to reach
 *             the start of the next limb in `res`.
 * @param mu The mean value of the distribution.
 * @param sigma The variance of the distribution.
 *  
 * @retval `-1` if an error occurs. In this case the error is from a syscall in rand_uniform.
 * @retval `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int normal_random_vec(int64_t limb_len, 
                      double* res, int64_t res_size, int64_t res_sl,
                      double mu, double sigma
){
    for(int i = 0; i < res_size; i++)
    {
        for(int j = 0; j < limb_len ; j++)
        {    
            if(rand_normal(res + (i * res_sl + j), mu, sigma) < 0) 
                return -1;
        }
    }

    return 0;
}
