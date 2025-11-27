#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <time.h>

int main(void)
{
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);

    gsl_rng_set(r, time(NULL));

    double x = gsl_rng_uniform(r);

    double mean = 0.0;
    double sigma = 1.0; // sigma = standard deviation
    double z = mean + gsl_ran_gaussian(r, sigma);  

    printf("Uniform  Number : %f\n", x);
    printf("Gaussian Number : %f\n", z);
    gsl_rng_free(r);
    return 0;
}