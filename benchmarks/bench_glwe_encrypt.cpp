#include <benchmark/benchmark.h>

extern "C" {
#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
}

#include "params.h"
#include "utils.hpp"

void bench_encrypt_rnx(benchmark::State& state)
{
	PvdaBackend* module = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	PolyUnivRnX* m                 = new_univ_rnx(params_glwe);
	GLWECiphertext* glwe_computed  = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, SKBITS);
	glwe_sk_prepare(module, sk_prep, sk);

	rnx_random_vec(m, params_glwe);

	for (auto _ : state)
	{
		glwe_secret_encrypt_rnx(module, glwe_computed, sk_prep, m);
		benchmark::DoNotOptimize(glwe_computed);
	}

	delete_univ_rnx(m);
	delete_glwe(glwe_computed);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
}

BENCHMARK(bench_encrypt_rnx);

void bench_encrypt_tnx(benchmark::State& state)
{
	PvdaBackend* module = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	PolyUnivTnX* m                 = new_univ_tnx(params_glwe);
	GLWECiphertext* glwe_computed  = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, SKBITS);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_random_pol_znx((PolyUniv*)m, NBASE, 64);

	for (auto _ : state)
	{
		glwe_secret_encrypt_tnx(module, glwe_computed, sk_prep, m);
		benchmark::DoNotOptimize(glwe_computed);
	}

	delete_univ_tnx(m);
	delete_glwe(glwe_computed);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
}

BENCHMARK(bench_encrypt_tnx);

void bench_encrypt_tnx_normalnoise(benchmark::State& state)
{
	PvdaBackend* module     = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_NORMAL);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	PolyUnivTnX* m                 = new_univ_tnx(params_glwe);
	GLWECiphertext* glwe_computed  = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, SKBITS);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_random_pol_znx((PolyUniv*)m, NBASE, 64);

	for (auto _ : state)
	{
		glwe_secret_encrypt_tnx(module, glwe_computed, sk_prep, m);
		benchmark::DoNotOptimize(glwe_computed);
	}

	delete_univ_tnx(m);
	delete_glwe(glwe_computed);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
}

// Benchmark ommited until normal noise implementation
//BENCHMARK(bench_encrypt_tnx_normalnoise);
