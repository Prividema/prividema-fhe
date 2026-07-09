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

void bench_enc_dec_rnx(benchmark::State& state)
{
	PvdaBackend* module = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	PolyUnivRnX* m                 = new_univ_rnx(params_glwe);
	GLWECiphertext* glwe_computed  = new_glwe(params_glwe);
	PolyBiv* result_biv            = new_biv(params_glwe);
	PolyUnivRnX* result_univ       = new_univ_rnx(params_glwe);

	uniform_glwe_secret_key(module, sk, SKBITS);
	glwe_sk_prepare(module, sk_prep, sk);

	rnx_random_vec(module, m, params_glwe);

	for (auto _ : state)
	{
		glwe_secret_encrypt_rnx(module, glwe_computed, sk_prep, m);
		glwe_secret_decrypt(module, result_biv, sk_prep, glwe_computed);
		biv_to_univ_rnx(params_glwe, result_univ, result_biv);
		benchmark::DoNotOptimize(result_univ);
	}

	delete_univ_rnx(m);
	delete_glwe(glwe_computed);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_biv(result_biv);
	delete_univ_rnx(result_univ);
}

BENCHMARK(bench_enc_dec_rnx);
