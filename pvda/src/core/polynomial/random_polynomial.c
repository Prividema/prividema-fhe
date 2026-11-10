#include "random_polynomial.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "backend_arithmetic.h"
#include "bivariate_polynomial.h"
#include "logger.h"
#include "univariate_polynomial.h"
#include "utils.h"

// =============================================
// |                                           |
// |                PolyUnivZnx                |
// |                                           |
// =============================================

int univ_znx_random_uniform_pow2(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n, uint64_t nb_bits)
{
	return pvda_rng_uniform_pow2_vec(module, res, n, nb_bits);
}

int univ_znx_random_uniform(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n, int64_t low, int64_t high)
{
	return pvda_rng_uniform_range_vec(module, res, n, low, high);
}

int univ_znx_random_binary(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n) { return pvda_rng_binary_vec(module, res, n); }

int univ_znx_random_ternary(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n)
{
	return pvda_rng_uniform_range_vec(module, res, n, -1, 1);
}

int univ_znx_random_binomial(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n, uint64_t eta)
{
	return pvda_rng_centered_binomial_vec(module, res, n, eta);
}

// Samples are drawn in chunks so that the temporary buffer stays small for large n.
#define NORMAL_CHUNK 1024

int univ_znx_random_normal(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n, double sigma)
{
	double tmp[NORMAL_CHUNK];
	for (uint64_t i = 0; i < n; i += NORMAL_CHUNK)
	{
		uint64_t len = n - i < NORMAL_CHUNK ? n - i : NORMAL_CHUNK;
		CHECK_CALL(pvda_rng_normal_vec(module, tmp, len, 0.0, sigma), "normal sampling failed in univ_znx_random_normal");
		for (uint64_t j = 0; j < len; ++j)
		{
			double r = round(tmp[j]);
			// [-2^63, 2^63) is exactly representable in double
			if (!(r >= -0x1p63 && r < 0x1p63)) RAISE_ERROR("Normal sample does not fit in an int64_t");
			res[i + j] = (int64_t)r;
		}
	}
	return 0;

cleanup:
	return -1;
}

// =============================================
// |                                           |
// |                PolyUnivRnX                |
// |                                           |
// =============================================

int univ_rnx_random_uniform(const PvdaBackend* module, PolyUnivRnX* res, uint64_t n, double low, double high)
{
	if (!isfinite(low) || !isfinite(high) || !(low < high))
		RAISE_ERROR("univ_rnx_random_uniform: low and high must be finite with low < high");

	// The 53 random bits of each sample are drawn as 64-bit words in res itself, then converted in place.
	CHECK_CALL(pvda_rng_bytes(module, res, n * sizeof(double)), "RNG failed in univ_rnx_random_uniform");
	for (uint64_t i = 0; i < n; ++i)
	{
		uint64_t w;
		memcpy(&w, res + i, sizeof(w));
		double u = (double)(w >> 11) * 0x1p-53;  // [0, 1)
		double v = low + (high - low) * u;
		res[i]   = v < high ? v : low;  // guard against rounding up to high
	}
	return 0;

cleanup:
	return -1;
}

int univ_rnx_random_normal(const PvdaBackend* module, PolyUnivRnX* res, uint64_t n, double mu, double sigma)
{
	return pvda_rng_normal_vec(module, res, n, mu, sigma);
}

// =============================================
// |                                           |
// |                PolyUnivTnX                |
// |                                           |
// =============================================

int univ_tnx_random_uniform(const PvdaBackend* module, PolyUnivTnX* res, uint64_t n)
{
	return pvda_rng_bytes(module, res, n * sizeof(PolyUnivTnX));
}

int univ_tnx_random_normal(const PvdaBackend* module, PolyUnivTnX* res, uint64_t n, double sigma)
{
	double tmp[NORMAL_CHUNK];
	for (uint64_t i = 0; i < n; i += NORMAL_CHUNK)
	{
		uint64_t len = n - i < NORMAL_CHUNK ? n - i : NORMAL_CHUNK;
		CHECK_CALL(pvda_rng_normal_vec(module, tmp, len, 0.0, sigma), "normal sampling failed in univ_tnx_random_normal");
		for (uint64_t j = 0; j < len; ++j)
		{
			// Reduce mod 1 into [-1/2, 1/2], then scale to 2^64 and round
			double r = tmp[j] - round(tmp[j]);
			double d = round(ldexp(r, 64));  // in [-2^63, 2^63]
			if (d >= 0x1p63) d -= 0x1p64;    // 1/2 == -1/2 mod 1
			res[i + j] = (uint64_t)(int64_t)d;
		}
	}
	return 0;

cleanup:
	return -1;
}

// =============================================
// |                                           |
// |               PolyUnivZnxMP               |
// |                                           |
// =============================================

int univ_znx_mp_random_uniform(const PvdaBackend* module, PolyUnivZnxMP* res)
{
	return pvda_rng_bytes(module, res->ptr, univ_znx_mp_bytes(res));
}

// =============================================
// |                                           |
// |                  PolyBiv                  |
// |                                           |
// =============================================

int biv_random_uniform_pow2(const PvdaBackend* module, PolyBiv* res, uint64_t nb_bits)
{
	for (uint64_t i = 0; i < res->nb_limbs; ++i)
		CHECK_CALL(pvda_rng_uniform_pow2_vec(module, res->ptr + i * res->stride, res->nn, nb_bits),
		           "uniform sampling failed in biv_random_uniform_pow2");
	return 0;

cleanup:
	return -1;
}

int biv_random_uniform(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res, uint64_t precision)
{
	if (precision > res->nb_limbs) RAISE_ERROR("biv_random_uniform: precision exceeds the number of limbs");

	PolyBiv head = biv_view(res->nn, precision, res->stride, res->ptr);
	return biv_random_uniform_pow2(module, &head, params_glwe->kappa);

cleanup:
	return -1;
}

int biv_random_normal(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res)
{
	int status = -1;

	PolyUnivRnX* rd_pol_univ = NULL;

	double sigma = glwe_params_stdev(params_glwe);
	if (isnan(sigma)) RAISE_ERROR("biv_random_normal: the noise type of the parameters has no standard deviation");

	rd_pol_univ = new_univ_rnx(params_glwe);
	CHECK_ALLOC(rd_pol_univ, "rd_pol_univ's malloc failed.");

	CHECK_CALL(univ_rnx_random_normal(module, rd_pol_univ, params_glwe->nn, 0.0, sigma),
	           "random normal vec generation failed");

	CHECK_CALL(univ_rnx_to_biv(params_glwe, res, rd_pol_univ, 0),
	           "univ to biv conversion failed in normal random biv poly generation");

	status = 0;

cleanup:
	delete_univ_rnx(rd_pol_univ);
	return status;
}

int biv_random_binomial(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res, uint64_t eta)
{
	uint64_t l_a = glwe_params_l_a(params_glwe);
	if (res->nb_limbs < l_a) RAISE_ERROR("biv_random_binomial: res has fewer than l_a limbs");

	for (uint64_t i = 0; i + 1 < l_a; ++i) memset(res->ptr + i * res->stride, 0, res->nn * sizeof(int64_t));
	return pvda_rng_centered_binomial_vec(module, res->ptr + (l_a - 1) * res->stride, res->nn, eta);

cleanup:
	return -1;
}

int biv_dft_random_uniform_pow2(const PvdaBackend* module, BE_PolyBivDFT* res_dft, uint64_t size,
                                uint64_t nb_bits)
{
	int status = -1;

	uint64_t nn        = pvda_module_extract_nn(module);
	int64_t* tmp_space = malloc(nn * size * sizeof(int64_t));
	CHECK_ALLOC(tmp_space, "malloc in biv_dft_random_uniform_pow2");

	PolyBiv tmp = biv_view(nn, size, (int64_t)nn, tmp_space);
	CHECK_CALL(biv_random_uniform_pow2(module, &tmp, nb_bits), "uniform sampling failed in biv_dft_random_uniform_pow2");
	CHECK_CALL(pvda_biv_to_dft(module, res_dft, size, &tmp), "DFT failed in biv_dft_random_uniform_pow2");

	status = 0;

cleanup:
	free(tmp_space);
	return status;
}

// =============================================
// |                                           |
// |                  Noise                    |
// |                                           |
// =============================================

/** @brief Samples the nn integer noise coefficients of the last limb, according to the noise type. */
static int sample_last_limb_noise(const PvdaBackend* module, const GLWEParams* params_glwe, PolyUnivZnx* err)
{
	uint64_t nn = params_glwe->nn;
	switch (params_glwe->noise_type)
	{
		case NOISE_UNIFORM_POWER_OF_TWO:
			if (params_glwe->fast_uniform_nb_bits == 0) break;
			return univ_znx_random_uniform_pow2(module, err, nn, params_glwe->fast_uniform_nb_bits);
		case NOISE_UNIFORM:
			return univ_znx_random_uniform(module, err, nn, -params_glwe->uniform_range, params_glwe->uniform_range);
		case NOISE_BINOMIAL:
			return univ_znx_random_binomial(module, err, nn, params_glwe->binomial_n / 2);
		case NOISE_NORMAL: {
			double rescaled = ldexp(params_glwe->normal_sigma, (int)params_glwe->kappa * (int)glwe_params_l_a(params_glwe));
			if (rescaled == 0) break;
			return univ_znx_random_normal(module, err, nn, rescaled);
		}
		default:
			return log_message(LOG_ERROR, "Noise type not implemented");
	}

	log_message(LOG_WARN, "Using error/noise with stdev 0");
	memset(err, 0, nn * sizeof(int64_t));
	return 0;
}

int biv_random_noise(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res)
{
	uint64_t l_a = glwe_params_l_a(params_glwe);
	if (res->nn != params_glwe->nn) RAISE_ERROR("biv_random_noise: res->nn differs from the parameters");
	if (res->nb_limbs < l_a) RAISE_ERROR("biv_random_noise: res has fewer than l_a limbs");

	for (uint64_t i = 0; i < res->nb_limbs; ++i)
		if (i != l_a - 1) memset(res->ptr + i * res->stride, 0, res->nn * sizeof(int64_t));
	return sample_last_limb_noise(module, params_glwe, res->ptr + (l_a - 1) * res->stride);

cleanup:
	return -1;
}
