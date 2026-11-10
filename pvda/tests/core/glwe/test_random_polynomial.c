#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "backend.h"
#include "backend_arithmetic.h"
#include "backend_rng.h"
#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "random_polynomial.h"
#include "univariate_polynomial.h"

// Tests of the library-side random polynomials, of the encryption noise types, and of the RNG of a backend being
// configurable independently of its arithmetic.

#define N 4096

/** @brief A backend whose RNG is seeded, so that the tests are reproducible. */
static PvdaBackend* new_test_module(uint64_t nn)
{
	PvdaSeed seed;
	for (int i = 0; i < PVDA_SEED_BYTES; ++i) seed.bytes[i] = (uint8_t)(3 * i + 1);
	PvdaBackend* module = pvda_new_spqlios_backend((int)nn);
	cr_assert_not_null(module);
	cr_assert(pvda_backend_set_rng(module, PVDA_RNG_CPU, &seed) == 0);
	return module;
}

static void mean_var(const double* x, uint64_t n, double* mean, double* var)
{
	double s = 0, s2 = 0;
	for (uint64_t i = 0; i < n; ++i) s += x[i];
	*mean = s / (double)n;
	for (uint64_t i = 0; i < n; ++i) s2 += (x[i] - *mean) * (x[i] - *mean);
	*var = s2 / (double)(n - 1);
}

static void mean_var_i64(const int64_t* x, uint64_t n, double* mean, double* var)
{
	double* d = malloc(n * sizeof(double));
	for (uint64_t i = 0; i < n; ++i) d[i] = (double)x[i];
	mean_var(d, n, mean, var);
	free(d);
}

// =============================================
// |                                           |
// |                PolyUnivZnx                |
// |                                           |
// =============================================

Test(random_polynomial, znx_small_supports)
{
	PvdaBackend* module = new_test_module(1024);
	int64_t x[N];
	int seen[3];

	cr_assert(univ_znx_random_binary(module, x, N) == 0);
	for (int i = 0; i < N; ++i) cr_assert(x[i] == 0 || x[i] == 1);

	memset(seen, 0, sizeof(seen));
	cr_assert(univ_znx_random_ternary(module, x, N) == 0);
	for (int i = 0; i < N; ++i)
	{
		cr_assert(x[i] >= -1 && x[i] <= 1);
		seen[x[i] + 1] = 1;
	}
	cr_assert(seen[0] && seen[1] && seen[2]);

	cr_assert(univ_znx_random_uniform(module, x, N, 10, 12) == 0);
	for (int i = 0; i < N; ++i) cr_assert(x[i] >= 10 && x[i] <= 12);

	cr_assert(univ_znx_random_uniform_pow2(module, x, N, 4) == 0);
	for (int i = 0; i < N; ++i) cr_assert(x[i] >= -8 && x[i] < 8);

	cr_assert(univ_znx_random_binomial(module, x, N, 3) == 0);
	for (int i = 0; i < N; ++i) cr_assert(x[i] >= -3 && x[i] <= 3);

	pvda_delete_backend(module);
}

Test(random_polynomial, znx_normal_is_rounded_gaussian)
{
	PvdaBackend* module = new_test_module(1024);
	int64_t* x   = malloc(16 * N * sizeof(int64_t));  // more than one internal chunk

	const double sigma = 40.0;
	cr_assert(univ_znx_random_normal(module, x, 16 * N, sigma) == 0);
	double mean, var;
	mean_var_i64(x, 16 * N, &mean, &var);
	cr_assert(fabs(mean) < 6 * sigma / sqrt(16.0 * N), "mean %g", mean);
	// Rounding adds a variance of about 1/12
	cr_assert(fabs(var / (sigma * sigma + 1.0 / 12) - 1) < 0.05, "var %g", var);

	// Samples that do not fit in an int64_t are rejected
	cr_assert(univ_znx_random_normal(module, x, 64, 1e30) == -1);
	cr_assert(univ_znx_random_normal(module, x, 64, -1.0) == -1);

	free(x);
	pvda_delete_backend(module);
}

// =============================================
// |                                           |
// |             PolyUnivRnX / TnX / MP        |
// |                                           |
// =============================================

Test(random_polynomial, rnx_uniform_and_normal)
{
	PvdaBackend* module = new_test_module(1024);
	double x[N];

	cr_assert(univ_rnx_random_uniform(module, x, N, -0.5, 0.5) == 0);
	for (int i = 0; i < N; ++i) cr_assert(x[i] >= -0.5 && x[i] < 0.5);
	double mean, var;
	mean_var(x, N, &mean, &var);
	cr_assert(fabs(var * 12 - 1) < 0.1, "var %g", var);

	cr_assert(univ_rnx_random_uniform(module, x, 4, 1.0, 1.0) == -1);
	cr_assert(univ_rnx_random_uniform(module, x, 4, 2.0, 1.0) == -1);
	cr_assert(univ_rnx_random_uniform(module, x, 4, 0.0, INFINITY) == -1);

	cr_assert(univ_rnx_random_normal(module, x, N, -2.0, 0.5) == 0);
	mean_var(x, N, &mean, &var);
	cr_assert(fabs(mean + 2.0) < 6 * 0.5 / sqrt((double)N), "mean %g", mean);
	cr_assert(fabs(var / 0.25 - 1) < 0.15, "var %g", var);

	pvda_delete_backend(module);
}

Test(random_polynomial, tnx_normal_wraps_on_the_torus)
{
	PvdaBackend* module = new_test_module(1024);
	PolyUnivTnX t[N];
	double x[N];

	// Small sigma: read back as signed numerators over 2^64
	const double sigma = ldexp(1.0, -20);
	cr_assert(univ_tnx_random_normal(module, t, N, sigma) == 0);
	for (int i = 0; i < N; ++i) x[i] = ldexp((double)(int64_t)t[i], -64);
	double mean, var;
	mean_var(x, N, &mean, &var);
	cr_assert(fabs(sqrt(var) / sigma - 1) < 0.1, "stdev %g", sqrt(var));

	// Huge sigma: still valid torus elements (no overflow), roughly uniform
	cr_assert(univ_tnx_random_normal(module, t, N, 1e6) == 0);
	int neg = 0;
	for (int i = 0; i < N; ++i) neg += (int64_t)t[i] < 0;
	cr_assert(neg > N / 2 - 300 && neg < N / 2 + 300);

	cr_assert(univ_tnx_random_uniform(module, t, N) == 0);
	neg = 0;
	for (int i = 0; i < N; ++i) neg += (int64_t)t[i] < 0;
	cr_assert(neg > N / 2 - 300 && neg < N / 2 + 300);

	pvda_delete_backend(module);
}

Test(random_polynomial, znx_mp_uniform_fills_every_word)
{
	PvdaBackend* module = new_test_module(1024);
	PolyUnivZnxMP* mp = new_univ_znx_mp(64, 5);
	cr_assert(univ_znx_mp_random_uniform(module, mp) == 0);

	// Every word position has its top bit set for some coefficient (probability of failure 2^-64 per word)
	for (uint64_t w = 0; w < 5; ++w)
	{
		int top = 0;
		for (uint64_t p = 0; p < 64; ++p) top |= (mp->ptr[p * 5 + w] >> 63) != 0;
		cr_assert(top, "word %lu", w);
	}

	delete_univ_znx_mp(mp);
	pvda_delete_backend(module);
}

// =============================================
// |                                           |
// |                  PolyBiv                  |
// |                                           |
// =============================================

Test(random_polynomial, biv_uniform_and_binomial)
{
	PvdaBackend* module = new_test_module(1024);
	GLWEParams* params_glwe = new_glwe_params(64, 1, 8, 4 * 2, 0, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t l_a            = glwe_params_l_a(params_glwe);
	PolyBiv* pol            = new_biv(params_glwe);

	// Limbs beyond precision are left untouched
	for (uint64_t i = 0; i < l_a * 64; ++i) pol->ptr[i] = 1000;
	cr_assert(biv_random_uniform(module, params_glwe, pol, 2) == 0);
	for (uint64_t i = 0; i < 2 * 64; ++i) cr_assert(pol->ptr[i] >= -128 && pol->ptr[i] < 128);
	for (uint64_t i = 2 * 64; i < l_a * 64; ++i) cr_assert(eq(i64, pol->ptr[i], 1000));
	cr_assert(biv_random_uniform(module, params_glwe, pol, l_a + 1) == -1);

	// Strided: only the nn coefficients of each limb are written
	int64_t buf[3 * 100];
	for (int i = 0; i < 300; ++i) buf[i] = 1000;
	PolyBiv strided = biv_view(64, 3, 100, buf);
	cr_assert(biv_random_uniform_pow2(module, &strided, 2) == 0);
	for (int l = 0; l < 3; ++l)
		for (int p = 0; p < 100; ++p)
		{
			int64_t v = buf[l * 100 + p];
			if (p < 64) cr_assert(v >= -2 && v < 2);
			else cr_assert(eq(i64, v, 1000));
		}

	// Binomial: last limb only
	cr_assert(biv_random_binomial(module, params_glwe, pol, 4) == 0);
	for (uint64_t i = 0; i < (l_a - 1) * 64; ++i) cr_assert(eq(i64, pol->ptr[i], 0));
	for (uint64_t i = (l_a - 1) * 64; i < l_a * 64; ++i) cr_assert(pol->ptr[i] >= -4 && pol->ptr[i] <= 4);

	delete_biv(pol);
	delete_glwe_params(params_glwe);
	pvda_delete_backend(module);
}

Test(random_polynomial, biv_normal_needs_a_standard_deviation)
{
	PvdaBackend* module = new_test_module(1024);

	// Power-of-2 uniform noise has no defined stdev
	GLWEParams* params_pow2 = new_glwe_params(1024, 1, 19, 3 * 2, ldexp(1.0, -40), NOISE_UNIFORM_POWER_OF_TWO);
	PolyBiv* pol            = new_biv(params_pow2);
	cr_assert(biv_random_normal(module, params_pow2, pol) == -1);
	delete_biv(pol);
	delete_glwe_params(params_pow2);

	const double sigma        = ldexp(1.0, -20);
	GLWEParams* params_normal = new_glwe_params(1024, 1, 19, 3 * 2, sigma, NOISE_NORMAL);
	pol                       = new_biv(params_normal);
	double x[1024];
	cr_assert(biv_random_normal(module, params_normal, pol) == 0);
	biv_to_univ_rnx(params_normal, x, pol, 0);
	double mean, var;
	mean_var(x, 1024, &mean, &var);
	cr_assert(fabs(sqrt(var) / sigma - 1) < 0.15, "stdev %g", sqrt(var));

	delete_biv(pol);
	delete_glwe_params(params_normal);
	pvda_delete_backend(module);
}

Test(random_polynomial, biv_dft_uniform)
{
	PvdaBackend* module = new_test_module(64);
	BE_PolyBivDFT* dft  = pvda_new_biv_dft(module, 3);
	int64_t back[3 * 64];
	PolyBiv back_biv = biv_view(64, 3, 64, back);

	cr_assert(biv_dft_random_uniform_pow2(module, dft, 3, 5) == 0);
	cr_assert(pvda_biv_idft(module, &back_biv, dft, 3) == 0);
	// The inverse DFT is unnormalized: coefficients are those of the sampled polynomials scaled by the backend.
	// Here we only check that something non-trivial was transformed.
	int nonzero = 0;
	for (int i = 0; i < 3 * 64; ++i) nonzero += back[i] != 0;
	cr_assert(nonzero > 0);
	cr_assert(biv_dft_random_uniform_pow2(module, dft, 3, 0) == -1);

	pvda_delete_biv_dft(module, dft);
	pvda_delete_backend(module);
}

// =============================================
// |                                           |
// |                  Noise                    |
// |                                           |
// =============================================

typedef struct
{
	NoiseType type;
	int check_variance;  // NOISE_UNIFORM_POWER_OF_TWO picks the smallest power of 2 covering sigma: no exact stdev
} NoiseCase;

static const NoiseCase noise_cases[] = {
    {NOISE_UNIFORM_POWER_OF_TWO, 0},
    {NOISE_UNIFORM, 1},
    {NOISE_BINOMIAL, 1},
    {NOISE_NORMAL, 1},
};

Test(noise, every_noise_type_on_last_limb)
{
	PvdaBackend* module = pvda_new_spqlios_backend(1024);

	// K * l_a = 57 bits, sigma = 2^-50: rescaled stdev 2^7 on the last limb
	for (size_t c = 0; c < sizeof(noise_cases) / sizeof(noise_cases[0]); ++c)
	{
		GLWEParams* params_glwe = new_glwe_params(1024, 1, 19, 3 * 2, ldexp(1.0, -50), noise_cases[c].type);
		cr_assert_not_null(params_glwe, "noise type %d", noise_cases[c].type);
		uint64_t l_a = glwe_params_l_a(params_glwe);

		PolyBiv* noisy = new_biv(params_glwe);
		uint64_t l     = noisy->nb_limbs;
		int64_t* last  = malloc(16 * 1024 * sizeof(int64_t));
		for (int r = 0; r < 16; ++r)
		{
			for (uint64_t i = 0; i < l * 1024; ++i) noisy->ptr[i] = 77;  // overwritten everywhere
			cr_assert(biv_random_noise(module, params_glwe, noisy) == 0, "noise type %d", noise_cases[c].type);
			for (uint64_t i = 0; i < (l_a - 1) * 1024; ++i) cr_assert(eq(i64, noisy->ptr[i], 0));
			for (uint64_t i = l_a * 1024; i < l * 1024; ++i) cr_assert(eq(i64, noisy->ptr[i], 0));
			memcpy(last + r * 1024, noisy->ptr + (l_a - 1) * 1024, 1024 * sizeof(int64_t));
		}

		double mean, var;
		mean_var_i64(last, 16 * 1024, &mean, &var);
		cr_assert(var > 0, "noise type %d", noise_cases[c].type);
		if (noise_cases[c].check_variance)
		{
			double expected = ldexp(glwe_params_stdev(params_glwe), 57);
			cr_assert(fabs(sqrt(var) / expected - 1) < 0.05, "noise type %d: stdev %g, expected %g",
			          noise_cases[c].type, sqrt(var), expected);
			cr_assert(fabs(sqrt(var) / 128.0 - 1) < 0.05, "noise type %d: stdev %g", noise_cases[c].type, sqrt(var));
		}

		free(last);
		delete_biv(noisy);
		delete_glwe_params(params_glwe);
	}
	pvda_delete_backend(module);
}

Test(noise, strided_and_invalid_outputs)
{
	PvdaBackend* module     = new_test_module(64);
	GLWEParams* params_glwe = new_glwe_params(64, 1, 19, 3 * 2, ldexp(1.0, -50), NOISE_NORMAL);
	uint64_t l_a            = glwe_params_l_a(params_glwe);

	// Strided output: only the nn coefficients of each limb are written
	int64_t buf[3 * 100];
	for (int i = 0; i < 300; ++i) buf[i] = 1000;
	PolyBiv strided = biv_view(64, l_a, 100, buf);
	cr_assert(biv_random_noise(module, params_glwe, &strided) == 0);
	int nonzero = 0;
	for (uint64_t l = 0; l < l_a; ++l)
		for (int p = 0; p < 100; ++p)
		{
			int64_t v = buf[l * 100 + p];
			if (p >= 64) cr_assert(eq(i64, v, 1000));
			else if (l + 1 < l_a) cr_assert(eq(i64, v, 0));
			else nonzero += v != 0;
		}
	cr_assert(nonzero > 0);

	// Too few limbs, wrong ring dimension
	PolyBiv short_biv = biv_view(64, l_a - 1, 64, buf);
	cr_assert(biv_random_noise(module, params_glwe, &short_biv) == -1);
	PolyBiv wrong_nn = biv_view(32, l_a, 32, buf);
	cr_assert(biv_random_noise(module, params_glwe, &wrong_nn) == -1);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(module);
}

Test(noise, binomial_rejects_huge_noise)
{
	// Rescaled stdev 2^17: n = 2^36 trials per sample is refused
	cr_assert_null(new_glwe_params(1024, 1, 19, 3 * 2, ldexp(1.0, -40), NOISE_BINOMIAL));
}

Test(noise, encrypt_decrypt_with_every_noise_type)
{
	PvdaBackend* module = pvda_new_spqlios_backend(1024);

	for (size_t c = 0; c < sizeof(noise_cases) / sizeof(noise_cases[0]); ++c)
	{
		const double sigma      = ldexp(1.0, -50);
		GLWEParams* params_glwe = new_glwe_params(1024, 1, 19, 3 * 2, sigma, noise_cases[c].type);

		GLWESecretKey* sk              = new_glwe_secret_key(params_glwe);
		GLWESecretKeyPrepared* sk_prep = new_glwe_secret_key_prepared(module, params_glwe);
		PolyBiv* m                     = new_biv(params_glwe);
		PolyBiv* phase                 = new_biv(params_glwe);
		PolyBiv* dec                   = new_biv(params_glwe);
		GLWECiphertext* ct             = new_glwe(params_glwe);
		double m_rnx[1024], dec_rnx[1024];

		cr_assert(glwe_sk_random_binary(module, sk) == 0);
		cr_assert(glwe_sk_prepare(module, sk_prep, sk) == 0);
		cr_assert(biv_random_uniform(module, params_glwe, m, 2) == 0);

		cr_assert(biv_random_noise(module, params_glwe, phase) == 0);
		cr_assert(biv_add(module, phase, phase, m) == 0);
		cr_assert(glwe_secret_encrypt_phase(module, ct, sk_prep, phase) == 0);
		cr_assert(glwe_secret_decrypt(module, dec, sk_prep, ct) == 0);

		biv_to_univ_rnx(params_glwe, m_rnx, m, 0);
		biv_to_univ_rnx(params_glwe, dec_rnx, dec, 0);
		// Every noise type has a stdev of about sigma (power-of-2 uniform: at most 2 sigma), so errors stay < 16 sigma
		for (int p = 0; p < 1024; ++p)
		{
			double err = fabs(remainder(dec_rnx[p] - m_rnx[p], 1.0));
			cr_assert(err < 16 * sigma, "noise type %d: coef %d error %g", noise_cases[c].type, p, err);
		}

		delete_glwe(ct);
		delete_biv(dec);
		delete_biv(phase);
		delete_biv(m);
		delete_glwe_secret_key_prepared(sk_prep);
		delete_glwe_secret_key(sk);
		delete_glwe_params(params_glwe);
	}
	pvda_delete_backend(module);
}

// =============================================
// |                                           |
// |        Seeded RNG of the backend          |
// |                                           |
// =============================================

/** @brief Generates a key and encrypts a fixed message with a backend whose RNG is seeded with the given type. */
static void seeded_encryption(PvdaRngType type, const PvdaSeed* seed, int64_t* sk_out, GLWECiphertext* ct_out,
                              GLWEParams* params_glwe)
{
	PvdaBackend* module = pvda_new_spqlios_backend((int)params_glwe->nn);
	cr_assert(pvda_backend_set_rng(module, type, seed) == 0);

	GLWESecretKey* sk              = new_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = new_glwe_secret_key_prepared(module, params_glwe);
	PolyBiv* phase                 = new_biv(params_glwe);

	cr_assert(glwe_sk_random_ternary(module, sk) == 0);
	cr_assert(glwe_sk_prepare(module, sk_prep, sk) == 0);
	cr_assert(biv_random_noise(module, params_glwe, phase) == 0);
	cr_assert(glwe_secret_encrypt_phase(module, ct_out, sk_prep, phase) == 0);
	memcpy(sk_out, sk->values, params_glwe->nn * params_glwe->k * sizeof(int64_t));

	delete_biv(phase);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_glwe_secret_key(sk);
	pvda_delete_backend(module);
}

Test(rng_in_module, same_seed_same_key_and_ciphertext)
{
	GLWEParams* params_glwe = new_glwe_params(256, 2, 19, 3 * 3, ldexp(1.0, -50), NOISE_NORMAL);
	uint64_t sk_len         = params_glwe->nn * params_glwe->k;

	for (int t = 0; t < 2; ++t)
	{
		PvdaRngType type = t ? PVDA_RNG_OPENSSL : PVDA_RNG_CPU;
		if (!pvda_rng_type_available(type)) continue;

		PvdaSeed seed, other;
		cr_assert(pvda_seed_generate(&seed) == 0);
		cr_assert(pvda_seed_generate(&other) == 0);

		int64_t* sk_a = malloc(sk_len * sizeof(int64_t));
		int64_t* sk_b = malloc(sk_len * sizeof(int64_t));
		int64_t* sk_c = malloc(sk_len * sizeof(int64_t));
		GLWECiphertext* ct_a = new_glwe(params_glwe);
		GLWECiphertext* ct_b = new_glwe(params_glwe);
		GLWECiphertext* ct_c = new_glwe(params_glwe);

		seeded_encryption(type, &seed, sk_a, ct_a, params_glwe);
		seeded_encryption(type, &seed, sk_b, ct_b, params_glwe);
		seeded_encryption(type, &other, sk_c, ct_c, params_glwe);

		cr_assert(memcmp(sk_a, sk_b, sk_len * sizeof(int64_t)) == 0, "type %d", type);
		cr_assert(memcmp(ct_a->vec, ct_b->vec, glwe_params_bytes(params_glwe)) == 0, "type %d", type);
		cr_assert(memcmp(ct_a->vec, ct_c->vec, glwe_params_bytes(params_glwe)) != 0, "type %d", type);

		delete_glwe(ct_a);
		delete_glwe(ct_b);
		delete_glwe(ct_c);
		free(sk_a);
		free(sk_b);
		free(sk_c);
	}
	delete_glwe_params(params_glwe);
}
