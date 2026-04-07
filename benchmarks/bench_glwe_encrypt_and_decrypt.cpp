#include <benchmark/benchmark.h>

#include <cmath>

extern "C" {
#include "bivariate_polynomial.h"
#include "glwe.h"
#include "glwe_transform_key.h"
#include "univariate_polynomial.h"
}

#define NBASE       (1 << 14)
#define KBASE       1
#define KAPPABASE   19
#define NLIMBSBASE  (15 * 2)
#define LBASE       NLIMBSBASE / (KBASE + 1)
#define SIGMABASE   -(LBASE / 2 + 1) * KAPPABASE

#define PROB_FACTOR 3

void test_benchmark(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	// Since the message is drawn in Zn[X,Y], there is no decomposition error. Thus, the error should be smaller than 3*sigma 99.73% of the time
	double err_length = 3 * sigma;

	//! Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKey* sk             = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft      = alloc_glwe_secret_key_dft(params_glwe);
	PolyUnivRnX* m                = new_univ_rnx(params_glwe);
	GLWECiphertext* glwe_computed = new_glwe(params_glwe);
	PolyBiv* result_biv           = new_biv_poly(params_glwe);
	PolyUnivRnX* result_univ      = new_univ_rnx(params_glwe);

	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	//The input message, for now sampled normally since we cannot sample uniformly in the torus right now
	normal_random_vec(m, NBASE, 0.0, 0.1);

	for (auto _ : state)
	{
		glwe_secret_encrypt_rnx(module, glwe_computed, sk_dft, m);
		glwe_secret_decrypt(module, result_biv, sk_dft, glwe_computed);
		biv_to_univ_rnx(params_glwe, result_univ, result_biv);
		benchmark::DoNotOptimize(result_univ);
	}

	delete_univ_rnx(m);
	delete_glwe(glwe_computed);
	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
	free(result_biv);
	delete_univ_rnx(result_univ);
}

BENCHMARK(test_benchmark);
