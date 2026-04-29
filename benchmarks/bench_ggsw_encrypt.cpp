#include <benchmark/benchmark.h>

#include <cmath>

extern "C" {
#include "bivariate_polynomial.h"
#include "ggsw_arithmetic.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
}

#define NBASE      (1 << 14)
#define KBASE      1
#define KAPPABASE  19
#define NLIMBSBASE (15 * 2)
#define LBASE      NLIMBSBASE / (KBASE + 1)

void test_ggsw_encrypt(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NLIMBSBASE);

	GLWESecretKey* sk             = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft      = alloc_glwe_secret_key_dft(params_glwe);
	PolyUniv* m                   = new_univ(params_glwe);
	GGSWCiphertext* ggsw_computed = new_ggsw(params_ggsw);
	PolyBiv* result_biv           = new_biv_poly(params_glwe);
	PolyUnivRnX* result_univ      = new_univ_rnx(params_glwe);

	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	uniform_random_vec(NBASE, m, 1, NBASE, 4);

	for (auto _ : state)
	{
		ggsw_secret_encrypt(module, ggsw_computed, sk_dft, m);
		benchmark::DoNotOptimize(ggsw_computed);
	}

	delete_univ(m);
	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
	free(result_biv);
	delete_univ_rnx(result_univ);
	delete_ggsw(ggsw_computed);
}

BENCHMARK(test_ggsw_encrypt);
