#include "rng.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <threads.h>

#include "logger.h"
#include "spqlios_alias.h"
#include "utils.h"

#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <bcrypt.h>
#include <windows.h>
#pragma comment(lib, "bcrypt.lib")
#endif

#ifdef __linux__
#include <sys/random.h>
#endif

// On some distros math.h doesn't define M_PI so we define it here just in case.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
    Read a random number depending on the OS :
    - On Windows : Uses Windows' Cryptographic API called CNG.
    - On MACOS/FreeBSD : Call to arc4random_buf.
                        According to arc4random's doc, the whole program crashes
                        if an error occurs during the generation.
    - On other Linux distributions : read /dev/urandom.
*/
int read_rand(uint64_t* result, size_t bytes)
{
// For Windows
#ifdef _WIN32
	NTSTATUS status = BCryptGenRandom(NULL, (PUCHAR)result, sizeof(*result), BCRYPT_USE_SYSTEM_PREFERRED_RNG);
	if (status != STATUS_SUCCESS) return log_message(LOG_ERROR, "BCryptGenRandom() Failed");

// For MACOS/FreeBSD
// According to arc4random's doc, the function crashes if an error occurs :
// "Cryptographic randomness is considered fundamental — if it’s broken, continuing execution is unsafe."
#elif defined(__APPLE__) || defined(__FreeBSD__)
	arc4random_buf(result, sizeof(*result));

// For other Linux Distro
#elif defined(__linux__)

	size_t rand_bytes = getrandom(result, bytes, 0);
	if (rand_bytes != bytes) return -1;
	return 0;

#else
	FILE* f = fopen("/dev/urandom", "rb");
	if (!f) return log_perror("fopen");
	int r = fread(result, sizeof(*result), 1, f);
	fclose(f);
	if (r != 1) return log_perror("fread");
#endif

	return 0;
}
static inline void reduce_uniform_n(int64_t* tgt, int n_bits)
{
	int shft = 64 - n_bits;
	*tgt     = ((*tgt) << shft) >> shft;
}

int rand_uniform(int64_t* result, uint64_t nb_bits)
{
	// As result points to an uint64_t  nb_bits shall not exceed its size
	assert(nb_bits <= 64);

	// If nb_bits equals the max. size, we just have to convert r to an int64_t.
	if (nb_bits == 8 * sizeof(int64_t))
		return read_rand((uint64_t*)result, 8);

	else
	{
		if (read_rand((uint64_t*)result, nb_bits / 8 + (nb_bits % 8 != 0)) < 0) return -1;

		reduce_uniform_n(result, nb_bits);

		return 1;
	}

	return 0;
}

/*
    Approximate inverse error function (erfinv) using Winitzki
    approximation

    If you have a uniform random variable X which is uniformly distributed
    between 0 and 1, you can map it to any distribution by using its inverse CDF.

    Since the iCDF of the normal distribution can be expressed with erfinv.
    We use this function for performance purpose.

    See Winitzki's paper called "A handy approximation for the error
    function and its inverse" or https://en.wikipedia.org/wiki/Error_function.
 */
double erfinv(double x)
{
	double a           = 0.147;
	double ln_1minusx2 = log(1.0 - x * x);
	double term1       = (2 / (M_PI * a)) + (ln_1minusx2 / 2.0);
	double term2       = ln_1minusx2 / a;
	return (x > 0 ? 1 : -1) * sqrt(sqrt(term1 * term1 - term2) - term1);
}

/*
    This function transforms a uniformly sampled variable into a normally
    distributed variable using the inverse Cumulative Distribution Function
    (CDF).
 */
int rand_normal(double* result, double mu, double sigma)
{
	// Generate a uniform number in [0, 2^64]
	uint64_t uniform;
	CHECK_CALL(read_rand(&uniform, 8), "Rng failed in rand_normal");

	// Scale uniform in (0,1) to U : U still follows a uniform distribution.
	double uu = ((double)uniform) / ((double)UINT64_MAX);

	// Compute Z the inverse CDF of the normal distribution applied to U.
	double zz = sqrt(2.0) * erfinv(2.0 * uu - 1.0);

	// Scale and Shift with mu and sigma.
	// Z follows a normal distribution in (0,1)
	// Thus result will follow (mu, sigma)
	*result = mu + sigma * zz;

	return 0;
cleanup:
	return -1;
}

int uniform_random_pol_znx(PolyUniv* res, uint64_t N, uint64_t nb_bits)
{
	CHECK_CALL(read_rand((uint64_t*)res, sizeof(int64_t) * N), "rng error");
	for (uint64_t p = 0; p < N; p++)
	{
		reduce_uniform_n(res + p, nb_bits);
	}
	return 0;
cleanup:
	return -1;
}

//TODO: fix int64_t number of limbs and stride length
int uniform_random_vec(uint64_t limb_len, int64_t* res, int64_t nb_limbs, int64_t res_sl, uint64_t nb_bits)
{
	for (uint64_t i = 0; i < nb_limbs; i++)
		CHECK_CALL(uniform_random_pol_znx(res + i * res_sl, limb_len, nb_bits), "unifrom random vec failed");
	return 0;
cleanup:
	return -1;
}

int uniform_random_vec_znx_dft(const MODULE* module, VecUnivDFT* result_dft, uint64_t vec_size, uint64_t nb_bits)
{
	int status = -1;

	// Variables
	int64_t* tmp_space = NULL;

	// The degree of the cyclotomic polynomial
	uint64_t nn = pvda_module_extract_nn(module);

	// Pointer to a uniformly drawn Zn[X] vector of size = vec_size
	tmp_space = malloc(nn * vec_size * sizeof(int64_t));
	CHECK_ALLOC(tmp_space, "malloc in new_uniform_random_vec_znx_dft");

	// Draws uniformly in Zn[X] the vector elements
	for (int i = 0; i < vec_size; i++)
		for (int p = 0; p < nn; p++)
			CHECK_CALL(rand_uniform(tmp_space + i * nn + p, nb_bits),
			           "rand_uniform failed in uniform_random_vec_znx_dft");

	// Computes the vector in the DFT domain
	pvda_vec_znx_dft(module, result_dft, vec_size, tmp_space, vec_size, nn);

	status = 0;

cleanup:
	free(tmp_space);

	return status;
}

int add_normal_random_vec(double* res, size_t vec_size, const double* vec, double mu, double sigma)
{
	// TODO: possible performance improvement here
	for (int i = 0; i < vec_size; i++)
	{
		double tmp;
		CHECK_CALL(rand_normal(&tmp, mu, sigma), "add_normal_random_vec failed");
		res[i] = vec[i] + tmp;
	}

	return 0;
cleanup:
	return -1;
}

//TODO: int64_t size
int normal_random_vec(double* res, int64_t res_size, double mu, double sigma)
{
	for (int i = 0; i < res_size; i++) CHECK_CALL(rand_normal(res + i, mu, sigma), "normal_random_vec failed");

	return 0;
cleanup:
	return -1;
}
