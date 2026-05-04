#include <benchmark/benchmark.h>

#include <cmath>

extern "C" {
#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
}

#define NBASE      (1 << 14)
#define KBASE      1
#define KAPPABASE  19
#define NLIMBSBASE (15 * 2)
#define LBASE      NLIMBSBASE / (KBASE + 1)

void test_encrypt_rnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, 4 - (LBASE)*KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_FAST_UNIFORM);

	GLWESecretKey* sk             = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft      = alloc_glwe_secret_key_dft(params_glwe);
	PolyUnivRnX* m                = new_univ_rnx(params_glwe);
	GLWECiphertext* glwe_computed = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	//The input message, for now sampled normally since we cannot sample uniformly in the torus right now
	normal_random_vec(m, NBASE, 0.0, 0.1);

	for (auto _ : state)
	{
		glwe_secret_encrypt_rnx(module, glwe_computed, sk_dft, m);
		benchmark::DoNotOptimize(glwe_computed);
	}

	delete_univ_rnx(m);
	delete_glwe(glwe_computed);
	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
}

BENCHMARK(test_encrypt_rnx);

void test_encrypt_tnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, 4 - (LBASE)*KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_FAST_UNIFORM);

	GLWESecretKey* sk             = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft      = alloc_glwe_secret_key_dft(params_glwe);
	PolyUnivTnX* m                = new_univ_tnx(params_glwe);
	GLWECiphertext* glwe_computed = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	//The input message, for now sampled normally since we cannot sample uniformly in the torus right now
	uniform_random_pol_znx((PolyUniv*)m, NBASE, 64);

	for (auto _ : state)
	{
		glwe_secret_encrypt_tnx(module, glwe_computed, sk_dft, m);
		benchmark::DoNotOptimize(glwe_computed);
	}

	delete_univ_tnx(m);
	delete_glwe(glwe_computed);
	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
}

BENCHMARK(test_encrypt_tnx);

void test_encrypt_tnx_normalnoise(benchmark::State& state)
{
	double sigma = ldexp(1.0, 4 - (LBASE)*KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_NORMAL);

	GLWESecretKey* sk             = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft      = alloc_glwe_secret_key_dft(params_glwe);
	PolyUnivTnX* m                = new_univ_tnx(params_glwe);
	GLWECiphertext* glwe_computed = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	//The input message, for now sampled normally since we cannot sample uniformly in the torus right now
	uniform_random_pol_znx((PolyUniv*)m, NBASE, 64);

	for (auto _ : state)
	{
		glwe_secret_encrypt_tnx(module, glwe_computed, sk_dft, m);
		benchmark::DoNotOptimize(glwe_computed);
	}

	delete_univ_tnx(m);
	delete_glwe(glwe_computed);
	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
}

BENCHMARK(test_encrypt_tnx_normalnoise);
