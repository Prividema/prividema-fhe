#include <benchmark/benchmark.h>

extern "C" {
#include "bivariate_polynomial.h"
#include "ggsw_arithmetic.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "glwegadget_ciphertext.h"
#include "rng.h"
#include "univariate_polynomial.h"
}

#include "params.h"
#define MSGBITS 4  //Bit size of the message coefficients

void bench_glwegad_encrypt(benchmark::State& state)
{
	PvdaBackend* module = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* params_glwegad = new_glwegadget_params(params_glwe, KAPPABASE, LBASE);

	GLWESecretKey* sk                      = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep         = alloc_glwe_secret_key_prepared(params_glwe);
	PolyUniv* m                            = new_univ(params_glwe);
	GLWEGadgetCiphertext* glwegad_computed = new_glwegadget(params_glwegad);
	PolyBiv* result_biv                    = new_biv(params_glwe);
	PolyUnivRnX* result_univ               = new_univ_rnx(params_glwe);

	uniform_glwe_secret_key(module, sk, SKBITS);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_pow2_random_pol_znx(m, NBASE, MSGBITS);

	for (auto _ : state)
	{
		glwegadget_secret_encrypt(module, glwegad_computed, sk_prep, m);
		benchmark::DoNotOptimize(glwegad_computed);
	}

	delete_univ(m);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_glwegadget_params(params_glwegad);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_biv(result_biv);
	delete_univ_rnx(result_univ);
	delete_glwegadget(glwegad_computed);
}

BENCHMARK(bench_glwegad_encrypt);
