#include "distributions.h"
//gsl_rng_uniform_int

int64_t* uniform_random_vec(int64_t* res, int64_t res_size);

int main(void)
{
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);

    gsl_rng_set(r, time(NULL));

    int x = gsl_rng_uniform_int(r,10);

    double mean = 0.0;
    double sigma = 1.0; // sigma = standard deviation
    double z = mean + gsl_ran_gaussian(r, sigma);  

    printf("Uniform  Number : %d\n", x);
    printf("Gaussian Number : %f\n", z);
    gsl_rng_free(r);
    return 0;
}