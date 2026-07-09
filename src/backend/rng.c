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

// On some distros math.h doesn't define M_PI so we define it here just in case.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
	return pvda_rand_uniform_pow2_vec(module, (int64_t*)res, nn, nb_bits);
}

int binary_random_pol_znx(const PvdaBackend* module, PolyUniv* res, uint64_t nn)
{
	return pvda_rand_uniform_binary_vec(module, (uint64_t*)res, nn);
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
