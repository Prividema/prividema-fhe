#include <benchmark/benchmark.h>

extern "C" {
#include "bivariate_polynomial.h"
#include "ggsw_arithmetic.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
}

#include "params.h"
#include "utils.hpp"
#define MSGBITS 12  //Bit size of m

void bench_unprepared_ggsw_ext_prod(benchmark::State& state)
{
	MODULE* module = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NLIMBSBASE);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	PolyUniv* m                    = new_univ(params_glwe);
	GLWECiphertext* glwe_input     = new_glwe(params_glwe);
	GLWECiphertext* glwe_computed  = new_glwe(params_glwe);
	GGSWCiphertext* ggsw           = new_ggsw(params_ggsw);
	PolyBiv* result_biv            = new_biv_poly(params_glwe);
	PolyUnivRnX* result_univ       = new_univ_rnx(params_glwe);
	PolyUnivRnX* m_glwe            = new_univ_rnx(params_glwe);

	uniform_glwe_secret_key(module, sk, SKBITS);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_random_vec(NBASE, m, 1, NBASE, MSGBITS);
	rnx_random_vec(m_glwe, params_glwe);
	ggsw_secret_encrypt(module, ggsw, sk_prep, m);
	glwe_secret_encrypt_rnx(module, glwe_input, sk_prep, m_glwe);

	for (auto _ : state)
	{
		ggsw_unprepared_external_product(module, glwe_computed, glwe_input, ggsw);
		benchmark::DoNotOptimize(glwe_computed);
	}

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_univ(m);
	delete_glwe(glwe_input);
	delete_glwe(glwe_computed);
	delete_ggsw(ggsw);
	free(result_biv);
	delete_univ_rnx(result_univ);
	delete_univ_rnx(m_glwe);

	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
}

BENCHMARK(bench_unprepared_ggsw_ext_prod);

void bench_ggsw_prepared_prod(benchmark::State& state)
{
	double sigma = ldexp(1.0, 4 - (LBASE)*KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NLIMBSBASE);

	GLWESecretKey* sk                 = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep    = alloc_glwe_secret_key_prepared(params_glwe);
	PolyUniv* m                       = new_univ(params_glwe);
	GLWECiphertext* glwe_input        = new_glwe(params_glwe);
	GLWECiphertext* glwe_computed     = new_glwe(params_glwe);
	GGSWCiphertext* ggsw              = new_ggsw(params_ggsw);
	GGSWCiphertextPrep* ggsw_prepared = new_ggsw_prep(params_ggsw);
	PolyBiv* result_biv               = new_biv_poly(params_glwe);
	PolyUnivRnX* result_univ          = new_univ_rnx(params_glwe);
	PolyUnivRnX* m_glwe               = new_univ_rnx(params_glwe);

	uniform_glwe_secret_key(module, sk, 3);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_random_vec(NBASE, m, 1, NBASE, 4);
	rnx_random_vec(m_glwe, params_glwe);
	ggsw_secret_encrypt(module, ggsw, sk_prep, m);
	ggsw_prepare(module, ggsw_prepared, ggsw);
	glwe_secret_encrypt_rnx(module, glwe_input, sk_prep, m_glwe);

	for (auto _ : state)
	{
		ggsw_external_product(module, glwe_computed, glwe_input, ggsw_prepared);
		benchmark::DoNotOptimize(glwe_computed);
	}

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_univ(m);
	delete_glwe(glwe_input);
	delete_glwe(glwe_computed);
	delete_ggsw(ggsw);
	delete_ggsw_prep(ggsw_prepared);
	free(result_biv);
	delete_univ_rnx(result_univ);
	delete_univ_rnx(m_glwe);

	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
}

BENCHMARK(bench_ggsw_prepared_prod);
