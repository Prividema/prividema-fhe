#include "distributions.h"
#include "utils.h"

/**
 * @brief Generates a Random Vector following an Uniform distribution.
 * 
 * @param res      The vector that stores the random data.
 * @param res_size The size of the Vector.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall in rand_uniform.
 * @retval - `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int uniform_random_vec(int64_t* res, int64_t res_size)
{
    for(int i = 0; i < res_size; i++)
        if(rand_uniform(res + i) < 0)
            return -1;
    return 0;
}