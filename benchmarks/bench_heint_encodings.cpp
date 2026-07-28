
#include <benchmark/benchmark.h>

#include <cstring>

extern "C" {
#include "backend.h"
#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "pvda_ffts.h"
#include "schemes/heint.h"
#include "univariate_polynomial.h"
}

#include "params.h"

void bench_heint_encoding(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivTnX* initial_vec = new_univ_tnx(params_glwe);
	PolyBiv* interm_vec      = new_biv(params_glwe);
	PolyUnivTnX* tmp_tnx     = new_univ_tnx(params_glwe);

	memset(initial_vec, 0, poly_univ_rnx_bytes(params_glwe));
	uint64_t nn = params_glwe->nn;

	generate_ntt_table(backend, TBASE);
	for (auto _ : state)
	{
		heint_encode(backend, params_glwe, interm_vec, nn, TBASE, initial_vec);
		benchmark::DoNotOptimize(interm_vec);
	}

	delete_univ_tnx(initial_vec);
	delete_biv(interm_vec);
	delete_univ_tnx(tmp_tnx);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(backend);
}
BENCHMARK(bench_heint_encoding);
/*

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

*/
