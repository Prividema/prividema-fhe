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


// Approximate inverse error function (erfinv) using Winitzki approximation
// Good enough for practical purposes
double erfinv(double x) {
    double a = 0.147;
    double ln_1minusx2 = log(1.0 - x*x);
    double term1 = (2 / (M_PI * a)) + (ln_1minusx2 / 2.0);
    double term2 = ln_1minusx2 / a;
    return (x > 0 ? 1 : -1) * sqrt( sqrt(term1 * term1 - term2) - term1 );
}

/**
 * Generates a random number following a normal distribution. mu = 0 , sigma = 
 * 
 * @param result A pointer that will point to the generated value.
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int rand_normal(int64_t *result, double sigma)
{
    // Generate an uniform number in [0, 2^64]
    uint64_t uniform;
    if (read_rand(&uniform) < 0)
        return -1;

    // Convert it in (0,1) with an uniform distribution
    double U = (uniform + 0.5) / ((double)UINT64_MAX);

    // Convert U in (0,1) with a normal distribution
    double Z = sqrt(2.0) * erfinv(2.0 * U - 1.0);

    // Convert Z in INT64 range
    double sigma = (double)INT64_MAX / 6.0;
    double Yd = sigma * Z;

    // Clamp to int64
    if (Yd < (double)INT64_MIN) Yd = (double)INT64_MIN;
    if (Yd > (double)INT64_MAX) Yd = (double)INT64_MAX;

    // Round to get the result
    *result = (int64_t)llround(Yd);

    return 0;
}