#include "ckks.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "backend.h"
#include "backend_private.h"  //todo remove
#include "logger.h"
#include "utils.h"

void encode_slow_internal(double* out, uint64_t n, double complex* in)
{
	for (uint64_t i = 0; i < n; ++i)
	{
		double complex sum = 0;
		uint64_t pw        = 1;
		uint64_t pw_mask   = 2 * n - 1;
		for (uint64_t j = 0; j < n / 2; ++j)
		{
			double complex expo = pw * i * (M_PI * I) / n;
			sum += cexp(-expo) * in[j];
			pw *= 5;
			pw &= pw_mask;
		}
		out[i] = 2.0 / n * creal(sum);
	}
}

void complex_bitrev(double complex* inout, uint64_t n)
{
	int64_t j = 0;
	for (int64_t i = 1; i < n; ++i)
	{
		uint64_t bit = n >> 1;

		while (j >= bit)
		{
			j -= bit;
			bit >>= 1;
		}

		j += bit;

		if (i < j)
		{
			double complex tmp = inout[i];
			inout[i]           = inout[j];
			inout[j]           = tmp;
		}
	}
}

int encode_ifft(const PvdaBackend* backend, double complex* inout, uint64_t slots)
{
	// Adapted from Lattigo's SpecialFFTDouble
	assert(slots >= 2);
	uint64_t k = slots;

	uint64_t log2k = next_pow2_log(k);
	uint64_t nn    = backend->pvda_fft_data->nn;
	uint64_t log2n = next_pow2_log(2 * nn);

	for (uint64_t lvl = log2k; lvl > 0; --lvl)
	{
		uint64_t lvl_len      = 1 << lvl;
		uint64_t half_lvl_len = lvl_len >> 1;
		uint64_t quad_lvl_len = lvl_len << 2;
		uint64_t log_gap      = log2n - 2 - lvl;
		uint64_t mask         = quad_lvl_len - 1;

		for (uint64_t i = 0; i < k; i += lvl_len)
		{
			for (uint64_t j = 0; j < half_lvl_len; ++j)
			{
				uint64_t p = i + j;
				double complex v1, v2;
				double complex root;

				v1 = inout[p] + inout[p + half_lvl_len];

				uint64_t root_idx = backend->pvda_fft_data->rotation_group[j] & mask;
				root_idx          = (quad_lvl_len - root_idx) << log_gap;
				root              = backend->pvda_fft_data->roots[root_idx];
				v2                = (inout[p] - inout[p + half_lvl_len]) * root;

				inout[p]                = v1;
				inout[p + half_lvl_len] = v2;
			}
		}
	}

	for (uint64_t i = 0; i < slots; ++i)
	{
		inout[i] /= (double)slots;
	}

	complex_bitrev(inout, slots);

	return 0;
}

int decode_ifft(const PvdaBackend* backend, double complex* inout, uint64_t slots)
{
	// Adapted from Lattigo's SpecialFFTDouble
	assert(slots >= 2);
	uint64_t k = slots;

	uint64_t log2k = next_pow2_log(k);
	uint64_t nn    = backend->pvda_fft_data->nn;
	uint64_t log2n = next_pow2_log(2 * nn);

	for (uint64_t lvl = log2k; lvl > 0; --lvl)
	{
		uint64_t lvl_len      = 1 << lvl;
		uint64_t half_lvl_len = lvl_len >> 1;
		uint64_t quad_lvl_len = lvl_len << 2;
		uint64_t log_gap      = log2n - 2 - lvl;
		uint64_t mask         = quad_lvl_len - 1;

		for (uint64_t i = 0; i < k; i += lvl_len)
		{
			for (uint64_t j = 0; j < half_lvl_len; ++j)
			{
				uint64_t p = i + j;
				double complex v1, v2;
				double complex root;

				v1 = inout[p] + inout[p + half_lvl_len];

				uint64_t root_idx = backend->pvda_fft_data->rotation_group[j] & mask;
				root_idx          = (quad_lvl_len - root_idx) << log_gap;
				root              = backend->pvda_fft_data->roots[root_idx];
				v2                = (inout[p] - inout[p + half_lvl_len]) * root;

				inout[p]                = v1;
				inout[p + half_lvl_len] = v2;
			}
		}
	}

	for (uint64_t i = 0; i < slots; ++i)
	{
		inout[i] /= (double)slots;
	}

	complex_bitrev(inout, slots);

	return 0;
}

int encode_internal(const PvdaBackend* backend, double* out, uint64_t slots, double complex* in, int inplace)
{
	int status = -1;
	if (!inplace)
	{
		double complex* tmp = calloc(slots, sizeof(double complex));
		if (!tmp)
		{
			log_message(LOG_ERROR, "Memory allocation failed in not-inplace CKKS encoding");
			return -1;
		}
		memcpy(tmp, in, slots * sizeof(double complex));
		in = tmp;
	}

	encode_ifft(backend, in, slots);

	for (int i = 0; i < slots; ++i)
	{
		out[i]         = creal(in[i]);
		out[i + slots] = cimag(in[i]);
	}
	status = 0;
cleanup:
	if (!inplace)
	{
		free(in);
	}
	return 0;
}

void decode_slow_internal(double complex* out, uint64_t n, double* in)
{
	uint64_t pw      = 1;
	uint64_t pw_mask = 2 * n - 1;
	for (uint64_t j = 0; j < n / 2; ++j)
	{
		double complex expo = pw * (M_PI * I) / n;
		double complex root = cexp(expo);
		double complex sum  = 0;
		for (int64_t i = (int64_t)n - 1; i >= 0; --i)
		{
			sum *= root;
			sum += in[i];
		}
		out[j] = sum;
		pw *= 5;
		pw &= pw_mask;
	}
}
