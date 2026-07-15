
#include <benchmark/benchmark.h>

#include <cstring>

extern "C" {
#include "backend.h"
#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "schemes/hefp.h"
#include "schemes/hefp_private.h"
#include "univariate_polynomial.h"
}

#include "params.h"

void bench_hefp_encoding(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyBiv* interm_vec      = new_biv(params_glwe);
	PolyUnivTnX* tmp_tnx     = new_univ_tnx(params_glwe);

	memset(initial_vec, 0, poly_univ_rnx_bytes(params_glwe));
	uint64_t nn = params_glwe->nn;

	for (auto _ : state)
	{
		hefp_encode(backend, params_glwe, interm_vec, nn / 2, 0, (_Complex double*)initial_vec);
		benchmark::DoNotOptimize(interm_vec);
	}

	delete_univ_rnx(initial_vec);
	delete_biv(interm_vec);
	delete_univ_tnx(tmp_tnx);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(backend);
}

BENCHMARK(bench_hefp_encoding);

void bench_hefp_decoding(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivRnX* final_vec = new_univ_rnx(params_glwe);
	PolyBiv* interm_vec    = new_biv(params_glwe);
	PolyUnivTnX* tmp_tnx   = new_univ_tnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	for (auto _ : state)
	{
		hefp_decode(backend, params_glwe, (_Complex double*)final_vec, nn / 2, 0, interm_vec);
		benchmark::DoNotOptimize(final_vec);
	}

	delete_univ_rnx(final_vec);
	delete_biv(interm_vec);
	delete_univ_tnx(tmp_tnx);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(backend);
}

BENCHMARK(bench_hefp_decoding);

void bench_hefp_decoding_fft(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivRnX* final_vec = new_univ_rnx(params_glwe);
	memset(final_vec, 0, poly_univ_rnx_bytes(params_glwe));

	uint64_t nn = params_glwe->nn;

	for (auto _ : state)
	{
		hefp_decode_fft(backend, (_Complex double*)final_vec, nn / 2);
		benchmark::DoNotOptimize(final_vec);
	}

	delete_univ_rnx(final_vec);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(backend);
}

BENCHMARK(bench_hefp_decoding_fft);

void bench_hefp_decoding_fft_slow(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	memset(initial_vec, 0, poly_univ_rnx_bytes(params_glwe));

	uint64_t nn = params_glwe->nn;

	for (auto _ : state)
	{
		hefp_decode_slow_internal((_Complex double*)final_vec, nn / 2, initial_vec);
		benchmark::DoNotOptimize(final_vec);
	}

	delete_univ_rnx(initial_vec);
	delete_univ_rnx(final_vec);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(backend);
}

BENCHMARK(bench_hefp_decoding_fft_slow);
