
#include <benchmark/benchmark.h>

#include <cmath>

extern "C" {
#include "bivariate_polynomial.h"
#include "ggsw_arithmetic.h"
#include "ggsw_params.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "glwegadget_arithmetic.h"
#include "rng.h"
#include "univariate_polynomial.h"
}

#include "params.h"

void test_glwegad_auto(benchmark::State& state)
{
	MODULE* module = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);
	GGSWParams* params_ggsw             = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NLIMBSBASE);
	GLWEGadgetParams* params_glwegadget = new_glwegadget_params(params_glwe, KAPPABASE, LBASE);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	GLWEAutomorphismKey* auto_key  = new_automorphism_key(params_glwegadget);

	PolyUnivTnX* m_univ_tnx     = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_expected_tnx = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_observed_tnx = new_univ_tnx(params_glwe);
	PolyBiv* m_auto             = new_biv(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);
	GLWECiphertext* glwe_res    = new_glwe(params_glwe);
	GLWECiphertext* glwe_norm   = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 3);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_random_pol_znx((PolyUniv*)m_univ_tnx, params_glwe->nn, 64);
	glwe_secret_encrypt_tnx(module, glwe_ct, sk_prep, m_univ_tnx);

	int auto_p = 7;

	compute_automorphism_key(module, auto_key, sk_prep, auto_p);

	for (auto _ : state)
	{
		glwegadget_automorphism(module, glwe_res, auto_key, glwe_ct);
		benchmark::DoNotOptimize(glwe_res);
	}

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_automorphism_key(auto_key);

	delete_univ_tnx(m_univ_tnx);
	delete_univ_tnx(m_observed_tnx);
	delete_univ_tnx(m_expected_tnx);
	delete_glwe(glwe_ct);
	delete_glwe(glwe_norm);
	delete_glwe(glwe_res);

	delete_biv(m_auto);

	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
	delete_glwegadget_params(params_glwegadget);
}

BENCHMARK(test_glwegad_auto);
