
#include <benchmark/benchmark.h>

#include <cmath>
#include <cstring>

extern "C" {
#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "glwegadget_arithmetic.h"
#include "glwegadget_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
}

#define NBASE      (1 << 12)
#define KBASE      1
#define KAPPABASE  19
#define NLIMBSBASE (4 * 2)
#define LBASE      NLIMBSBASE / (KBASE + 1)

void test_expand_trance(benchmark::State& state)
{
	double sigma = ldexp(1.0, 4 - (LBASE)*KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_FAST_UNIFORM);
	GLWEGadgetParams* params_glwegadget = new_glwegadget_params(params_glwe, KAPPABASE, LBASE);
	params_glwe->fast_uniform_nb_bits   = 0;
	sigma                               = 0;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);

	PolyUnivRnX* m_univ_rnx     = new_univ_rnx(params_glwe);
	PolyUnivTnX* m_expected_tnx = new_univ_tnx(params_glwe);
	PolyUnivRnX* m_observed_rnx = new_univ_rnx(params_glwe);
	PolyBiv* biv_tmp            = new_biv_poly(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	GLWEAutomorphismKSK** ksks = (GLWEAutomorphismKSK**)calloc(2 * params_glwe->nn, sizeof(GLWEAutomorphismKSK*));
	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		ksks[p]   = new_automorphism_ksk(params_glwegadget);
		prepare_automorphism_key(module, ksks[p], sk_prep, p);
	}
	int bund = NBASE;

	memset(m_univ_rnx, 0, poly_univ_bytes(params_glwe));

	// Get the message in univariate RnX form for expected result
	normal_random_vec(m_univ_rnx, bund, 0, 0.01);

	glwe_secret_encrypt_rnx(module, glwe_ct, sk_prep, m_univ_rnx);

	GLWECiphertext** results = (GLWECiphertext**)calloc(bund, sizeof(GLWECiphertext*));
	for (int i = 0; i < bund; ++i)
	{
		results[i] = new_glwe(params_glwe);
	}

	for (auto _ : state)
	{
		glwegadget_trace_expand(module, results, bund, glwe_ct, (const GLWEAutomorphismKSK**)ksks, 2 * params_glwe->nn);

		benchmark::DoNotOptimize(results);
	}

	int a = 1;
	for (int i = 0; i < bund; ++i)
	{
		delete_glwe(results[i]);
	}
	free(results);

	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		delete_automorphism_ksk(ksks[p]);
	}
	free(ksks);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	free(biv_tmp);
	delete_univ_rnx(m_univ_rnx);
	delete_univ_rnx(m_observed_rnx);
	delete_univ_tnx(m_expected_tnx);
	delete_glwe(glwe_ct);

	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_glwegadget_params(params_glwegadget);
}

BENCHMARK(test_expand_trance);
