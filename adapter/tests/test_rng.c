#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "backend.h"
#include "backend_rng.h"
#include "impls/rng_cpu.h"
#include "stat_utils.h"

// Tests of the RNG of a backend: configuration, seeding, reference vectors and every distribution on every RNG type.
//
// Statistical checks use tolerances of about 6 standard deviations, so that they fail with negligible probability
// on a correct implementation.

#define NB_SAMPLES 100000

// =============================================
// |                                           |
// |              Helpers                      |
// |                                           |
// =============================================

/** @brief Every RNG configuration available in this build: each backend, in system and seeded mode. */
typedef struct
{
	PvdaRngType type;
	int seeded;
} RngConfig;

static const RngConfig all_configs[] = {
    {PVDA_RNG_CPU, 0},
    {PVDA_RNG_CPU, 1},
    {PVDA_RNG_OPENSSL, 0},
    {PVDA_RNG_OPENSSL, 1},
};
#define NB_CONFIGS (sizeof(all_configs) / sizeof(all_configs[0]))

static PvdaSeed fixed_seed(uint8_t start)
{
	PvdaSeed seed;
	for (int i = 0; i < PVDA_SEED_BYTES; ++i) seed.bytes[i] = (uint8_t)(start + i);
	return seed;
}

/** @brief Creates a backend whose RNG is (type, seed), or returns NULL if the RNG type is not compiled in. */
static PvdaBackend* new_rng_module(PvdaRngType type, const PvdaSeed* seed)
{
	if (!pvda_rng_type_available(type)) return NULL;
	PvdaBackend* module = pvda_new_spqlios_backend(16);
	cr_assert_not_null(module);
	cr_assert(pvda_backend_set_rng(module, type, seed) == 0, "RNG setup failed for type %d", type);
	return module;
}

/** @brief Creates a backend with the RNG of a configuration, or returns NULL if its type is not compiled in. */
static PvdaBackend* new_config_module(RngConfig c)
{
	PvdaSeed seed = fixed_seed(7);
	return new_rng_module(c.type, c.seeded ? &seed : NULL);
}

#define FOR_EACH_RNG(rng)                                                   \
	for (size_t cfg_i = 0; cfg_i < NB_CONFIGS; ++cfg_i)                     \
		for (PvdaBackend* rng = new_config_module(all_configs[cfg_i]); rng; \
		     pvda_delete_backend(rng), rng = NULL)

static void mean_var_i64(const int64_t* x, uint64_t n, double* mean, double* var)
{
	double s = 0, s2 = 0;
	for (uint64_t i = 0; i < n; ++i) s += (double)x[i];
	*mean = s / (double)n;
	for (uint64_t i = 0; i < n; ++i) s2 += ((double)x[i] - *mean) * ((double)x[i] - *mean);
	*var = s2 / (double)(n - 1);
}

// =============================================
// |                                           |
// |              Creation                     |
// |                                           |
// =============================================

Test(rng_backend, cpu_always_available) { cr_assert(eq(int, pvda_rng_type_available(PVDA_RNG_CPU), 1)); }

Test(rng_backend, set_rng_matches_availability)
{
	PvdaSeed seed       = fixed_seed(0);
	PvdaBackend* module = pvda_new_spqlios_backend(16);
	for (size_t i = 0; i < NB_CONFIGS; ++i)
	{
		RngConfig c         = all_configs[i];
		PvdaRngType before  = pvda_backend_rng_type(module);
		int before_seeded   = pvda_backend_rng_is_seeded(module);
		int status          = pvda_backend_set_rng(module, c.type, c.seeded ? &seed : NULL);
		if (pvda_rng_type_available(c.type))
		{
			cr_assert(eq(int, status, 0));
			cr_assert(eq(int, pvda_backend_rng_type(module), c.type));
			cr_assert(eq(int, pvda_backend_rng_is_seeded(module), c.seeded));
		}
		else
		{
			// The previous RNG is kept
			cr_assert(eq(int, status, -1));
			cr_assert(eq(int, pvda_backend_rng_type(module), before));
			cr_assert(eq(int, pvda_backend_rng_is_seeded(module), before_seeded));
		}
		uint8_t b[8];
		cr_assert(pvda_rng_bytes(module, b, sizeof(b)) == 0);
	}
	pvda_delete_backend(module);
}

Test(rng_backend, default_rng_of_a_new_backend)
{
	PvdaRngType expected = pvda_rng_type_available(PVDA_RNG_OPENSSL) ? PVDA_RNG_OPENSSL : PVDA_RNG_CPU;
	cr_assert(eq(int, pvda_default_rng_type(), expected));

	PvdaBackend* module = pvda_new_spqlios_backend(16);
	cr_assert(eq(int, pvda_backend_rng_type(module), expected));
	cr_assert(eq(int, pvda_backend_rng_is_seeded(module), 0));
	uint8_t b[8];
	cr_assert(pvda_rng_bytes(module, b, sizeof(b)) == 0);
	pvda_delete_backend(module);
}

Test(rng_backend, seed_generate_gives_distinct_seeds)
{
	PvdaSeed a, b;
	cr_assert(pvda_seed_generate(&a) == 0);
	cr_assert(pvda_seed_generate(&b) == 0);
	cr_assert(memcmp(a.bytes, b.bytes, PVDA_SEED_BYTES) != 0);
}

Test(rng_backend, reseeding_restarts_the_stream)
{
	PvdaSeed seed       = fixed_seed(9);
	PvdaBackend* module = new_rng_module(PVDA_RNG_CPU, &seed);
	uint8_t first[100], again[100];
	cr_assert(pvda_rng_bytes(module, first, sizeof(first)) == 0);
	cr_assert(pvda_backend_set_rng(module, PVDA_RNG_CPU, &seed) == 0);
	cr_assert(pvda_rng_bytes(module, again, sizeof(again)) == 0);
	cr_assert(memcmp(first, again, sizeof(first)) == 0);

	// Back to system mode
	cr_assert(pvda_backend_set_rng(module, PVDA_RNG_CPU, NULL) == 0);
	cr_assert(eq(int, pvda_backend_rng_is_seeded(module), 0));
	pvda_delete_backend(module);
}

// =============================================
// |                                           |
// |          Seeded mode: reference vectors   |
// |                                           |
// =============================================

// RFC 8439, appendix A.1, test vector #1: all-zero key and nonce, block counter 0.
static const uint8_t chacha20_zero_block0[64] = {
    0x76, 0xb8, 0xe0, 0xad, 0xa0, 0xf1, 0x3d, 0x90, 0x40, 0x5d, 0x6a, 0xe5, 0x53, 0x86, 0xbd, 0x28,
    0xbd, 0xd2, 0x19, 0xb8, 0xa0, 0x8d, 0xed, 0x1a, 0xa8, 0x36, 0xef, 0xcc, 0x8b, 0x77, 0x0d, 0xc7,
    0xda, 0x41, 0x59, 0x7c, 0x51, 0x57, 0x48, 0x8d, 0x77, 0x24, 0xe0, 0x3f, 0xb8, 0xd8, 0x4a, 0x37,
    0x6a, 0x43, 0xb8, 0xf4, 0x15, 0x18, 0xa1, 0x1c, 0xc3, 0x87, 0xb6, 0x69, 0xb2, 0xee, 0x65, 0x86};

Test(rng_cpu, chacha20_reference_vector)
{
	const uint32_t key[8] = {0};
	uint8_t out[64];
	pvda_chacha20_block(key, 0, 0, out);
	cr_assert(memcmp(out, chacha20_zero_block0, 64) == 0);
}

Test(rng_cpu, seeded_stream_is_chacha20_keystream)
{
	PvdaSeed seed;
	memset(seed.bytes, 0, PVDA_SEED_BYTES);
	PvdaBackend* rng = new_rng_module(PVDA_RNG_CPU, &seed);

	uint8_t out[128], block1[64];
	cr_assert(pvda_rng_bytes(rng, out, sizeof(out)) == 0);
	cr_assert(memcmp(out, chacha20_zero_block0, 64) == 0);

	const uint32_t key[8] = {0};
	pvda_chacha20_block(key, 1, 0, block1);
	cr_assert(memcmp(out + 64, block1, 64) == 0);

	pvda_delete_backend(rng);
}

Test(rng_cpu, seed_is_little_endian_key)
{
	PvdaSeed seed = fixed_seed(0);  // bytes 0, 1, ..., 31
	PvdaBackend* rng = new_rng_module(PVDA_RNG_CPU, &seed);

	uint32_t key[8];
	for (int i = 0; i < 8; ++i)
		key[i] = (uint32_t)(4 * i) | (uint32_t)(4 * i + 1) << 8 | (uint32_t)(4 * i + 2) << 16 | (uint32_t)(4 * i + 3) << 24;
	uint8_t expected[64], out[64];
	pvda_chacha20_block(key, 0, 0, expected);
	cr_assert(pvda_rng_bytes(rng, out, sizeof(out)) == 0);
	cr_assert(memcmp(out, expected, 64) == 0);

	pvda_delete_backend(rng);
}

Test(rng_openssl, seeded_stream_is_aes256_ctr_keystream)
{
	if (!pvda_rng_type_available(PVDA_RNG_OPENSSL)) cr_skip_test("OpenSSL backend not compiled in");

	// AES-256 with the all-zero key on the all-zero block (FIPS-197 / NIST known answer)
	static const uint8_t expected[16] = {0xdc, 0x95, 0xc0, 0x78, 0xa2, 0x40, 0x89, 0x89,
	                                     0xad, 0x48, 0xa2, 0x14, 0x92, 0x84, 0x20, 0x87};
	PvdaSeed seed;
	memset(seed.bytes, 0, PVDA_SEED_BYTES);
	PvdaBackend* rng = new_rng_module(PVDA_RNG_OPENSSL, &seed);

	uint8_t out[16];
	cr_assert(pvda_rng_bytes(rng, out, sizeof(out)) == 0);
	cr_assert(memcmp(out, expected, 16) == 0);

	pvda_delete_backend(rng);
}

// =============================================
// |                                           |
// |          Seeded mode: determinism         |
// |                                           |
// =============================================

Test(rng_seeded, same_seed_same_stream_whatever_the_split)
{
	for (size_t t = 0; t < 2; ++t)
	{
		PvdaRngType type = t ? PVDA_RNG_OPENSSL : PVDA_RNG_CPU;
		if (!pvda_rng_type_available(type)) continue;

		PvdaSeed seed = fixed_seed(42);
		PvdaBackend* a = new_rng_module(type, &seed);
		PvdaBackend* b = new_rng_module(type, &seed);

		uint8_t out_a[1000], out_b[1000];
		cr_assert(pvda_rng_bytes(a, out_a, sizeof(out_a)) == 0);
		// Same stream when requested in odd-sized pieces
		cr_assert(pvda_rng_bytes(b, out_b, 1) == 0);
		cr_assert(pvda_rng_bytes(b, out_b + 1, 63) == 0);
		cr_assert(pvda_rng_bytes(b, out_b + 64, 100) == 0);
		cr_assert(pvda_rng_bytes(b, out_b + 164, 836) == 0);
		cr_assert(memcmp(out_a, out_b, sizeof(out_a)) == 0, "type %d", type);

		pvda_delete_backend(a);
		pvda_delete_backend(b);
	}
}

Test(rng_seeded, same_seed_same_distribution_samples)
{
	for (size_t t = 0; t < 2; ++t)
	{
		PvdaRngType type = t ? PVDA_RNG_OPENSSL : PVDA_RNG_CPU;
		if (!pvda_rng_type_available(type)) continue;

		PvdaSeed seed = fixed_seed(3);
		PvdaBackend* a = new_rng_module(type, &seed);
		PvdaBackend* b = new_rng_module(type, &seed);

		int64_t xa[257], xb[257];
		double na[33], nb[33];
		cr_assert(pvda_rng_uniform_range_vec(a, xa, 257, -5, 1000) == 0);
		cr_assert(pvda_rng_uniform_range_vec(b, xb, 257, -5, 1000) == 0);
		cr_assert(memcmp(xa, xb, sizeof(xa)) == 0);
		cr_assert(pvda_rng_normal_vec(a, na, 33, 1.0, 2.0) == 0);
		cr_assert(pvda_rng_normal_vec(b, nb, 33, 1.0, 2.0) == 0);
		cr_assert(memcmp(na, nb, sizeof(na)) == 0);

		pvda_delete_backend(a);
		pvda_delete_backend(b);
	}
}

Test(rng_seeded, different_seeds_and_backends_differ)
{
	PvdaSeed s1 = fixed_seed(1), s2 = fixed_seed(2);
	PvdaBackend* a = new_rng_module(PVDA_RNG_CPU, &s1);
	PvdaBackend* b = new_rng_module(PVDA_RNG_CPU, &s2);
	uint8_t out_a[64], out_b[64];
	cr_assert(pvda_rng_bytes(a, out_a, 64) == 0);
	cr_assert(pvda_rng_bytes(b, out_b, 64) == 0);
	cr_assert(memcmp(out_a, out_b, 64) != 0);
	pvda_delete_backend(a);
	pvda_delete_backend(b);

	if (pvda_rng_type_available(PVDA_RNG_OPENSSL))
	{
		a = new_rng_module(PVDA_RNG_CPU, &s1);
		b = new_rng_module(PVDA_RNG_OPENSSL, &s1);
		cr_assert(pvda_rng_bytes(a, out_a, 64) == 0);
		cr_assert(pvda_rng_bytes(b, out_b, 64) == 0);
		cr_assert(memcmp(out_a, out_b, 64) != 0);
		pvda_delete_backend(a);
		pvda_delete_backend(b);
	}
}

// =============================================
// |                                           |
// |              Distributions                |
// |                                           |
// =============================================

Test(rng_distributions, bytes_cover_all_values)
{
	FOR_EACH_RNG(rng)
	{
		uint8_t buf[8192];
		int seen[256] = {0};
		cr_assert(pvda_rng_bytes(rng, buf, sizeof(buf)) == 0);
		for (size_t i = 0; i < sizeof(buf); ++i) seen[buf[i]] = 1;
		for (int v = 0; v < 256; ++v) cr_assert(seen[v], "byte %d never drawn", v);
		cr_assert(pvda_rng_bytes(rng, buf, 0) == 0);
	}
}

Test(rng_distributions, uniform_pow2_range_and_coverage)
{
	int64_t* x                = malloc(NB_SAMPLES * sizeof(int64_t));
	static const uint64_t b[] = {1, 3, 17, 63};
	FOR_EACH_RNG(rng)
	{
		for (size_t k = 0; k < sizeof(b) / sizeof(b[0]); ++k)
		{
			int64_t lo = -(INT64_C(1) << (b[k] - 1)), hi = (INT64_C(1) << (b[k] - 1)) - 1;
			cr_assert(pvda_rng_uniform_pow2_vec(rng, x, NB_SAMPLES, b[k]) == 0);
			int64_t mn = INT64_MAX, mx = INT64_MIN;
			for (uint64_t i = 0; i < NB_SAMPLES; ++i)
			{
				cr_assert(x[i] >= lo && x[i] <= hi, "nb_bits %lu: %ld out of range", b[k], x[i]);
				mn = x[i] < mn ? x[i] : mn;
				mx = x[i] > mx ? x[i] : mx;
			}
			// Small ranges: both bounds must be reached
			if (b[k] <= 3)
			{
				cr_assert(eq(i64, mn, lo));
				cr_assert(eq(i64, mx, hi));
			}
			double mean, var;
			mean_var_i64(x, NB_SAMPLES, &mean, &var);
			double sd = ldexp(1.0, (int)b[k] - 1) / sqrt(3.0);  // stdev of the uniform
			cr_assert(fabs(mean + 0.5) < 6 * sd / sqrt(NB_SAMPLES) + 1, "nb_bits %lu: mean %g", b[k], mean);
		}
		// 64 bits: both signs appear
		cr_assert(pvda_rng_uniform_pow2_vec(rng, x, 1000, 64) == 0);
		int neg = 0;
		for (int i = 0; i < 1000; ++i) neg += x[i] < 0;
		cr_assert(neg > 400 && neg < 600);
	}
	free(x);
}

Test(rng_distributions, uniform_pow2_rejects_invalid_bits)
{
	int64_t x[4];
	FOR_EACH_RNG(rng)
	{
		cr_assert(pvda_rng_uniform_pow2_vec(rng, x, 4, 0) == -1);
		cr_assert(pvda_rng_uniform_pow2_vec(rng, x, 4, 65) == -1);
	}
}

Test(rng_distributions, uniform_range)
{
	int64_t* x = malloc(NB_SAMPLES * sizeof(int64_t));
	FOR_EACH_RNG(rng)
	{
		// Small range: every value is drawn, frequencies are balanced (chi-squared, 7 degrees of freedom)
		uint64_t counts[8] = {0};
		cr_assert(pvda_rng_uniform_range_vec(rng, x, NB_SAMPLES, -3, 4) == 0);
		for (uint64_t i = 0; i < NB_SAMPLES; ++i)
		{
			cr_assert(x[i] >= -3 && x[i] <= 4);
			counts[x[i] + 3]++;
		}
		double chi2 = 0, expected = NB_SAMPLES / 8.0;
		for (int v = 0; v < 8; ++v) chi2 += (counts[v] - expected) * (counts[v] - expected) / expected;
		cr_assert(chi2 < 30.0, "chi2 = %g", chi2);  // P(chi2_7 > 30) < 1e-4

		// Non power-of-2 width, close to INT64 limits
		cr_assert(pvda_rng_uniform_range_vec(rng, x, 1000, INT64_MAX - 4, INT64_MAX) == 0);
		for (int i = 0; i < 1000; ++i) cr_assert(x[i] >= INT64_MAX - 4);
		cr_assert(pvda_rng_uniform_range_vec(rng, x, 1000, INT64_MIN, INT64_MIN + 2) == 0);
		for (int i = 0; i < 1000; ++i) cr_assert(x[i] <= INT64_MIN + 2);

		// Degenerate range, full range, invalid range
		cr_assert(pvda_rng_uniform_range_vec(rng, x, 100, 17, 17) == 0);
		for (int i = 0; i < 100; ++i) cr_assert(eq(i64, x[i], 17));
		cr_assert(pvda_rng_uniform_range_vec(rng, x, 100, INT64_MIN, INT64_MAX) == 0);
		cr_assert(pvda_rng_uniform_range_vec(rng, x, 100, 1, 0) == -1);
	}
	free(x);
}

Test(rng_distributions, binary)
{
	int64_t* x = malloc(NB_SAMPLES * sizeof(int64_t));
	FOR_EACH_RNG(rng)
	{
		// Odd size to exercise the last partial word
		uint64_t n = NB_SAMPLES - 3;
		cr_assert(pvda_rng_binary_vec(rng, x, n) == 0);
		double mean, var;
		for (uint64_t i = 0; i < n; ++i) cr_assert(x[i] == 0 || x[i] == 1);
		mean_var_i64(x, n, &mean, &var);
		cr_assert(fabs(mean - 0.5) < 6 * 0.5 / sqrt((double)n), "mean %g", mean);
	}
	free(x);
}

Test(rng_distributions, centered_binomial)
{
	int64_t* x                  = malloc(NB_SAMPLES * sizeof(int64_t));
	static const uint64_t etas[] = {1, 2, 5, 64, 70};  // 70: more than one word per half-sample
	FOR_EACH_RNG(rng)
	{
		for (size_t k = 0; k < sizeof(etas) / sizeof(etas[0]); ++k)
		{
			int64_t eta = (int64_t)etas[k];
			cr_assert(pvda_rng_centered_binomial_vec(rng, x, NB_SAMPLES, etas[k]) == 0);
			for (uint64_t i = 0; i < NB_SAMPLES; ++i) cr_assert(x[i] >= -eta && x[i] <= eta);

			double mean, var, expected_var = eta / 2.0;
			mean_var_i64(x, NB_SAMPLES, &mean, &var);
			cr_assert(fabs(mean) < 6 * sqrt(expected_var / NB_SAMPLES), "eta %ld: mean %g", eta, mean);
			// The variance estimator has relative stdev ~ sqrt(2 / n) ~ 0.45%
			cr_assert(fabs(var / expected_var - 1) < 0.03, "eta %ld: var %g", eta, var);
		}
		// eta = 0: only zeros
		cr_assert(pvda_rng_centered_binomial_vec(rng, x, 100, 0) == 0);
		for (int i = 0; i < 100; ++i) cr_assert(eq(i64, x[i], 0));
	}
	free(x);
}

Test(rng_distributions, normal)
{
	double* x = malloc(NB_SAMPLES * sizeof(double));
	FOR_EACH_RNG(rng)
	{
		const double mu = 3.0, sigma = 0.25;
		// Odd size: the last Box-Muller pair is only half used
		uint64_t n = NB_SAMPLES - 1;
		cr_assert(pvda_rng_normal_vec(rng, x, n, mu, sigma) == 0);

		double s = 0, s2 = 0;
		for (uint64_t i = 0; i < n; ++i) s += x[i];
		double mean = s / n;
		for (uint64_t i = 0; i < n; ++i) s2 += (x[i] - mean) * (x[i] - mean);
		double var = s2 / (n - 1);

		cr_assert(fabs(mean - mu) < 6 * sigma / sqrt((double)n), "mean %g", mean);
		cr_assert(fabs(var / (sigma * sigma) - 1) < 0.03, "var %g", var);
		// Normality: Jarque-Bera statistic follows chi2_2, P(JB > 25) < 4e-6
		cr_assert(jarque_bera(x, (int)n) < 25.0, "Jarque-Bera %g", jarque_bera(x, (int)n));

		// sigma = 0 gives mu exactly
		cr_assert(pvda_rng_normal_vec(rng, x, 10, mu, 0.0) == 0);
		for (int i = 0; i < 10; ++i) cr_assert(eq(dbl, x[i], mu));
	}
	free(x);
}

Test(rng_distributions, normal_rejects_invalid_parameters)
{
	double x[4];
	FOR_EACH_RNG(rng)
	{
		cr_assert(pvda_rng_normal_vec(rng, x, 4, 0.0, -1.0) == -1);
		cr_assert(pvda_rng_normal_vec(rng, x, 4, 0.0, NAN) == -1);
		cr_assert(pvda_rng_normal_vec(rng, x, 4, INFINITY, 1.0) == -1);
	}
}
