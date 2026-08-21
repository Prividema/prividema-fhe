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
#define MSGBITS 4  //Bit size of the message coefficients

void bench_ggsw_encrypt(benchmark::State& state)
{
	PvdaBackend* module = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NLIMBSBASE);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	PolyUniv* m                    = new_univ(params_glwe);
	GGSWCiphertext* ggsw_computed  = new_ggsw(params_ggsw);
	PolyBiv* result_biv            = new_biv(params_glwe);
	PolyUnivRnX* result_univ       = new_univ_rnx(params_glwe);

	uniform_glwe_secret_key(module, sk, SKBITS);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_pow2_random_vec(module, NBASE, m, 1, NBASE, MSGBITS);

	for (auto _ : state)
	{
		ggsw_secret_encrypt(module, ggsw_computed, sk_prep, m);
		benchmark::DoNotOptimize(ggsw_computed);
	}

	delete_univ(m);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_biv(result_biv);
	delete_univ_rnx(result_univ);
	delete_ggsw(ggsw_computed);
}

BENCHMARK(bench_ggsw_encrypt);
