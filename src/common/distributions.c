#include "distributions.h"

/**
 * @brief Generates a Random Vector following an Uniform distribution.
 * 
 * @param res      The vector that stores the random data.
 * @param res_size The size of the Vector.
 * @param r        The RNG.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall in rand_seed.
 * @retval - `0` otherwise.
 * 
 * @note For each call to this function, we'll generate a new seed.
 */
int uniform_random_vec(int64_t* res, int64_t res_size, gsl_rng *r)
{
    if(gsl_refresh_seed(r) < 0)
        return -1;

    for(int i = 0; i < res_size; i++)
        res[i] = (int64_t) gsl_rng_get(r);
    return 0;
}