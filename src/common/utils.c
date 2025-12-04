#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M_PI 3.14159265358979323846

/**
 * Read a random number on /dev/urandom
 * 
 * @param result A pointer that will point to the generated value.
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int read_rand(uint64_t *result)
{
    FILE *f = fopen("/dev/urandom", "rb");
    if (!f) {
        perror("fopen");
        return -1;
    }

    if (fread(result, sizeof(result), 1, f) != 1) {
        perror("fread");
        fclose(f);
        return -1;
    }

    fclose(f);
    return 0;
}

/**
 * Generates a random number following an uniform distribution.
 * 
 * @param result A pointer that will point to the generated value.
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int rand_uniform(int64_t *result) { return (int64_t)read_rand(result); }


/**
 * @brief Approximate inverse error function (erfinv) using Winitzki approximation
 * 
 * @param x A uniformly random value between 0 and 1.
 * 
 * @note If you have a uniform random variable X which is uniformly distributed between 0 and 1, you can map it to any distribution by using its inverse CDF.
 * @note Since the iCDF of the normal distribution can be expressed with erfinv. We use this function for performance purpose.
 * @note See Winitzki's paper called "A handy approximation for the error function and its inverse" or https://en.wikipedia.org/wiki/Error_function.
 */
double erfinv(double x) {
    double a = 0.147;
    double ln_1minusx2 = log(1.0 - x*x);
    double term1 = (2 / (M_PI * a)) + (ln_1minusx2 / 2.0);
    double term2 = ln_1minusx2 / a;
    return (x > 0 ? 1 : -1) * sqrt( sqrt(term1 * term1 - term2) - term1 );
}

/**
 * Generates a random number between 0 and 1 that follows a normal distribution with given mu and sigma.
 * 
 * @param result A pointer that will point to the generated value.
 * @param mu     The mean value.
 * @param sigma  The variance.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 * 
 * @note This function transforms a uniformly sampled variable into a normally distributed variable using the inverse Cumulative Distribution Function (CDF).
 */
int rand_normal(double *result, double mu, double sigma)
{
    // Generate an uniform number in [0, 2^64]
    uint64_t uniform;
    if (read_rand(&uniform) < 0)
        return -1;

    // Scale uniform in (0,1) to U : U still follows an uniform distribution.
    double U = (uniform + 0.5) / ((double)UINT64_MAX);

    // Compute Z the inverse CDF of the normal distribution applied to U.
    double Z = sqrt(2.0) * erfinv(2.0 * U - 1.0);

    // Scale and Shift with mu and sigma.
    // result follow a normal distribution in (0,1)
    double result = mu + sigma * Z;

    return 0;
}