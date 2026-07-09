#include "rng.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <threads.h>

#include "backend.h"
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
inline int cpu_read_rand(uint64_t* result, size_t bytes)
{
// For Windows
#ifdef _WIN32
	NTSTATUS status = BCryptGenRandom(NULL, (PUCHAR)result, bytes, BCRYPT_USE_SYSTEM_PREFERRED_RNG);
	if (status != STATUS_SUCCESS) return log_message(LOG_ERROR, "BCryptGenRandom() Failed");

// For MACOS/FreeBSD
// According to arc4random's doc, the function crashes if an error occurs :
// "Cryptographic randomness is considered fundamental — if it’s broken, continuing execution is unsafe."
#elif defined(__APPLE__) || defined(__FreeBSD__)
	arc4random_buf(result, bytes);

// For Linux
#elif defined(__linux__)

	size_t rand_bytes = getrandom(result, bytes, 0);
	if (rand_bytes != bytes) return -1;
	return 0;

// I don't know what system this block below would be (and it would be quite slow)
// But for now we leave it for compatibility
#else
	// THIS IS VERY SLOW!
	FILE* f = fopen("/dev/urandom", "rb");
	if (!f) return log_perror("fopen");
	int r = fread(result, 1, bytes, f);
	fclose(f);
	if (r != bytes) return log_perror("fread");
#endif

	return 0;
}
static inline void reduce_uniform_n(int64_t* tgt, int n_bits)
{
	int shft = 64 - n_bits;
	*tgt     = (int64_t)((uint64_t)(*tgt) << shft) >> shft;
}

int ref_rand_uniform_pow2(const PvdaBackend* module, int64_t* result, uint64_t nb_bits)
{
	// As result points to an uint64_t  nb_bits shall not exceed its size
	assert(nb_bits <= 8 * sizeof(int64_t));

	// If nb_bits equals the max. size, we just have to convert r to an int64_t.
	if (nb_bits == 8 * sizeof(int64_t))
		return cpu_read_rand((uint64_t*)result, 8);

	else
	{
		if (cpu_read_rand((uint64_t*)result, INT_ROUND_UP_DIV(nb_bits, 8)) < 0) return -1;

		reduce_uniform_n(result, (int)nb_bits);

		return 1;
	}

	return 0;
}

int ref_rand_uniform(const PvdaBackend* module, int64_t* result, int64_t limit_down, int64_t limit_up)
{
	uint64_t max_delta = (uint64_t)limit_up - (uint64_t)limit_down;
	if (max_delta == UINT64_MAX) return -1;
	uint64_t bits = next_pow2_log(max_delta + 1);
	uint64_t mask = bits == 64 ? (UINT64_MAX) : (1ull << bits) - 1;

	uint64_t tmp = 0;

	int st;
	do
	{
		st = cpu_read_rand(&tmp, INT_ROUND_UP_DIV(bits, 8));
		if (st < 0) return -1;
		tmp &= mask;
	} while (tmp > max_delta);
	*result = (int64_t)((uint64_t)limit_down + tmp);
	return 0;
}

int rand_normal(double* result, double mu, double sigma)
{
	RAISE_ERROR("Normal number generation not implemented");
	/*
	// This function used to transforms a uniformly sampled variable into a normally
	// distributed variable using the inverse Cumulative Distribution Function
	// (CDF).
	// Its implementation has been removed since it used an approximation of the iCDF
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
  */
cleanup:
	return -1;
}

int uniform_pow2_random_pol_znx(const PvdaBackend* module, PolyUniv* res, uint64_t nn, uint64_t nb_bits)
{
	CHECK_CALL(cpu_read_rand((uint64_t*)res, sizeof(int64_t) * nn), "rng error");
	for (uint64_t p = 0; p < nn; p++)
	{
		reduce_uniform_n(res + p, (int)nb_bits);
	}
	return 0;
cleanup:
	return -1;
}

int binary_random_pol_znx(const PvdaBackend* module, PolyUniv* res, uint64_t nn)
{
	CHECK_CALL(cpu_read_rand((uint64_t*)res, sizeof(int64_t) * nn), "rng error");
	for (uint64_t p = 0; p < nn; p++)
	{
		res[p] &= 1;
	}
	return 0;
cleanup:
	return -1;
}

int uniform_random_pol_znx(const PvdaBackend* module, PolyUniv* res, uint64_t nn, int64_t low_bound, int64_t high_bound)
{
	for (uint64_t p = 0; p < nn; p++)
	{
		CHECK_CALL(pvda_rand_uniform(module, &res[p], low_bound, high_bound), "uniform RNG failed");
	}
	return 0;
cleanup:
	return -1;
}

int uniform_pow2_random_vec(const PvdaBackend* module, uint64_t limb_len, int64_t* res, uint64_t nb_limbs,
                            uint64_t res_sl, uint64_t nb_bits)
{
	for (uint64_t i = 0; i < nb_limbs; i++)
		CHECK_CALL(uniform_pow2_random_pol_znx(module, res + i * res_sl, limb_len, nb_bits),
		           "uniform random vec failed");
	return 0;
cleanup:
	return -1;
}

//Forward declaration due to legacy folder structure
PolyBiv new_biv_view(uint64_t nn, uint64_t l, int64_t stride, PolyBivUnderlying* ptr);

int uniform_random_vec_znx_dft(const PvdaBackend* module, VecUnivDFT* result_dft, uint64_t vec_size, uint64_t nb_bits)
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
			CHECK_CALL(pvda_rand_uniform_pow2(module, tmp_space + i * nn + p, nb_bits),
			           "rand_uniform failed in uniform_random_vec_znx_dft");

	// Computes the vector in the DFT domain
	PolyBiv tmp_biv = new_biv_view(nn, vec_size, nn, tmp_space);
	pvda_vec_znx_dft(module, result_dft, vec_size, &tmp_biv);

	status = 0;

cleanup:
	free(tmp_space);

	return status;
}

int add_normal_random_vec(const PvdaBackend* module, double* res, size_t vec_size, const double* vec, double mu,
                          double sigma)
{
	// Possible performance improvement here
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

int normal_random_vec(const PvdaBackend* module, double* res, uint64_t res_size, double mu, double sigma)
{
	for (int i = 0; i < res_size; i++) CHECK_CALL(rand_normal(res + i, mu, sigma), "normal_random_vec failed");

	return 0;
cleanup:
	return -1;
}
