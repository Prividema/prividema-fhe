#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Initializes an gsl Random Number Generator.
 * 
 * @param r The RNG.
 */
void gsl_init(gsl_rng *r)
{
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    r = gsl_rng_alloc(T);
}

/**
 * @brief Generates a random seed.
 * 
 * @param seed A pointer to the seed that has to be filled.
 *
 * @retval - `-1` if an error occurs. In this case, the error comes from the syscalls, and perror is called.
 * @retval - `0` otherwise.
 * 
 * @note This function is for now internal as we don't expose it in the header file. If needed, put it in the header.
 */
int rand_seed(unsigned long *seed)
{
    FILE *f = fopen("/dev/urandom", "rb");

    if (!f) {
        perror("fopen");
        return -1;
    }

    if (fread(seed, sizeof(seed), 1, f) != 1) {
        perror("fread");
        fclose(f);
        return -1;
    }
    fclose(f);
    return 0;
}

/**
 * @brief Refreshes the seed of the RNG.
 * 
 * @param r The RNG.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall in rand_seed.
 * @retval - `0` otherwise.
 * 
 * @note The seed is generated randomly by reading `/dev/urandom`.
 */
int gsl_refresh_seed(gsl_rng *r)
{
    unsigned long seed = 0;
    if(rand_seed(&seed) < 0)
        return -1;
    
    gsl_rng_set(r, seed);
    return 0;
}