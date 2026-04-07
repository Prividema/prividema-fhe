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
#define NLIMBSBASE  (20 * 2)
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
	GLWESecretKey* sk             = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft      = alloc_glwe_secret_key_dft(NBASE, KBASE);
	PolyBiv* m                    = new_biv_poly(params_glwe);
	PolyBiv* err                  = new_biv_poly(params_glwe);
	PolyBiv* phase                = new_biv_poly(params_glwe);
	GLWECiphertext* glwe_computed = new_glwe(params_glwe);

	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// The input message uniformly drawn in Zn[X,Y]
	uniform_random_biv_poly(params_glwe, m, LBASE / 2);

	// The input error normaly drawn in Zn[X,Y]
	normal_random_biv_poly(params_glwe, err);

	//! Computes with functions
	// Computes the message in Tn[X] with the base-2Kappa

	// The final phase = m + err
	add_biv_poly(params_glwe, phase, NBASE, m, NBASE, err, NBASE);

	for (auto _ : state)
	{
		glwe_secret_masking(module, glwe_computed, sk_dft, phase);
		benchmark::DoNotOptimize(glwe_computed);
	}

	free(m);
	free(err);
	free(phase);
	delete_glwe(glwe_computed);
	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
}

BENCHMARK(test_benchmark);
