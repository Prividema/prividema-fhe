#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
    #define _USE_MATH_DEFINES
    #include <bcrypt.h>
    #include <windows.h>
    #pragma comment(lib, "bcrypt.lib")
#endif

// On some distros math.h doesn't define M_PI so we define it here just in case.
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/**
 * Read a random number depending on the OS :
 * - On Windows : Uses Windows' Cryptographic API called CNG.
 * - On MACOS/FreeBSD : Call to arc4random_buf. 
 *                      According to arc4random's doc, the whole program crashes 
 *                      if an error occurs during the generation.
 * - On other Linux distributions : read /dev/urandom.
 */
int read_rand(uint64_t* result) {
    // For Windows
    #ifdef _WIN32
        NTSTATUS status = BCryptGenRandom(NULL, (PUCHAR)result, sizeof(*result), BCRYPT_USE_SYSTEM_PREFERRED_RNG);
        if (status != STATUS_SUCCESS) {
            fprintf(stderr, "Random generation failed\n");
            return -1;
        }

    // For MACOS/FreeBSD
    // According to arc4random's doc, the function crashes if an error occurs :
    // "Cryptographic randomness is considered fundamental — if it’s broken, continuing execution is unsafe."
    #elif defined(__APPLE__) || defined(__FreeBSD__)
        arc4random_buf(result, sizeof(*result));

    // For other Linux Distro
    #else
        FILE* f = fopen("/dev/urandom", "rb");
        if (!f) {
            perror("fopen");
            return -1;
        }

        if (fread(result, sizeof(*result), 1, f) != 1) {
            perror("fread");
            fclose(f);
            return -1;
        }

        fclose(f);
    #endif

    return 0;
}

int rand_uniform(int64_t* result, int nb_bits) {
    // As result points to an int64_t nb_bits shall not exceed its size
    if(nb_bits > 8 * sizeof(int64_t)) {
        fprintf(stderr, "rand_uniform() : nb_bits exceeds the maximum value (%ld).\n", 8 * sizeof(int64_t));
        return -1;
    }

    // Generate a random uint64_t
    // r is in the interval [0, UINT64_MAX]
    uint64_t r;
    int res;
    if((res = read_rand(&r)) < 0) 
        return -1;

    // If nb_bits equals the max. size, we just have to convert r to an int64_t.
    if (nb_bits == 8 * sizeof(int64_t))
        *result = (int64_t)r;

    // If nb_bits is not the max. size, r is in the interval [0, UINT64_MAX]
    // We bring r into the inteval [0, 2^nb_bits) with a modulo 
    // that is equivalent to truncating bits so we keep the cryptosafe property.
    // Then we apply an offset to get a result in [-2^(nb_bits-1), 2^(nb_bits-1))
    else {
        // Reduce modulo p = 2^nb_bits with a mask (1 << nb_bits) - 1
        // As r is still an unsigned int, it is now in [0, p) 
        uint64_t p = (1 << nb_bits);
        r &= p - 1;

        // Apply an offset so result is in [-p/2, p/2)
        *result = (int64_t)r - (int64_t)p/2;
    }

    return 0;
}

/**
 * Approximate inverse error function (erfinv) using Winitzki
 * approximation
 *
 * If you have a uniform random variable X which is uniformly distributed
 * between 0 and 1, you can map it to any distribution by using its inverse CDF.
 * 
 * Since the iCDF of the normal distribution can be expressed with erfinv.
 * We use this function for performance purpose.
 * 
 * See Winitzki's paper called "A handy approximation for the error
 * function and its inverse" or https://en.wikipedia.org/wiki/Error_function.
 */
double erfinv(double x) {
    double a = 0.147;
    double ln_1minusx2 = log(1.0 - x * x);
    double term1 = (2 / (M_PI * a)) + (ln_1minusx2 / 2.0);
    double term2 = ln_1minusx2 / a;
    return (x > 0 ? 1 : -1) * sqrt(sqrt(term1 * term1 - term2) - term1);
}

/**
 * This function transforms a uniformly sampled variable into a normally
 * distributed variable using the inverse Cumulative Distribution Function
 * (CDF).
 */
int rand_normal(double* result, double mu, double sigma) {
    // Generate an uniform number in [0, 2^64]
    uint64_t uniform;
    if (read_rand(&uniform) < 0) return -1;

    // Scale uniform in (0,1) to U : U still follows an uniform distribution.
    double U = (uniform + 0.5) / ((double)UINT64_MAX);

    // Compute Z the inverse CDF of the normal distribution applied to U.
    double Z = sqrt(2.0) * erfinv(2.0 * U - 1.0);

    // Scale and Shift with mu and sigma.
    // result follow a normal distribution in (0,1)
    *result = mu + sigma * Z;
    
    return 0;
}
