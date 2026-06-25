#include <benchmark/benchmark.h>

extern "C" {
#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "rng.h"
#include "univariate_polynomial.h"
}
#include "params.h"
#include "utils.hpp"

void bench_univ_biv_rnx(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivRnX* m = new_univ_rnx(params_glwe);
	PolyBiv* m_biv = new_biv(params_glwe);
	rnx_random_vec(m, params_glwe);

	for (auto _ : state)
	{
		univ_rnx_to_biv(params_glwe, m_biv, m, 0);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_rnx(m);
	delete_biv(m_biv);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(backend);
}

BENCHMARK(bench_univ_biv_rnx);

void bench_biv_normalize(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivRnX* m = new_univ_rnx(params_glwe);
	PolyBiv* m_biv = new_biv(params_glwe);

	rnx_random_vec(m, params_glwe);

	univ_rnx_to_biv(params_glwe, m_biv, m, 0);

	for (auto _ : state)
	{
		pvda_vec_znx_normalize_base2k(backend, params_glwe->kappa, m_biv, m_biv);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_rnx(m);
	delete_biv(m_biv);

	pvda_delete_backend(backend);
	delete_glwe_params(params_glwe);
}

BENCHMARK(bench_biv_normalize);

void bench_univ_biv_rnx_via_tnx(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivRnX* m     = new_univ_rnx(params_glwe);
	PolyBiv* m_biv     = new_biv(params_glwe);
	PolyUnivTnX* m_tnx = new_univ_tnx(params_glwe);
	rnx_random_vec(m, params_glwe);

	for (auto _ : state)
	{
		univ_rnx_to_tnx(params_glwe, m_tnx, m);
		univ_tnx_to_biv(params_glwe, m_biv, m_tnx, 0);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_rnx(m);
	delete_univ_tnx(m_tnx);
	delete_biv(m_biv);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(backend);
}

BENCHMARK(bench_univ_biv_rnx_via_tnx);

void bench_univ_biv_tnx(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivTnX* m = new_univ_tnx(params_glwe);
	PolyBiv* m_biv = new_biv(params_glwe);
	uniform_random_pol_znx((PolyUniv*)m, NBASE, 64);

	for (auto _ : state)
	{
		univ_tnx_to_biv(params_glwe, m_biv, m, 0);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_tnx(m);
	delete_biv(m_biv);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(backend);
}

BENCHMARK(bench_univ_biv_tnx);

void bench_biv_univ_tnx(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivTnX* m      = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_back = new_univ_tnx(params_glwe);
	PolyBiv* m_biv      = new_biv(params_glwe);
	uniform_random_pol_znx((PolyUniv*)m, NBASE, 64);
	univ_tnx_to_biv(params_glwe, m_biv, m, 0);

	for (auto _ : state)
	{
		biv_to_univ_tnx(params_glwe, m_back, m_biv);
		benchmark::DoNotOptimize(m_back);
	}

	delete_univ_tnx(m_back);
	delete_univ_tnx(m);
	delete_biv(m_biv);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(backend);
}

BENCHMARK(bench_biv_univ_tnx);

void bench_biv_univ_rnx(benchmark::State& state)
{
	PvdaBackend* backend = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);

	PolyUnivRnX* m      = new_univ_rnx(params_glwe);
	PolyUnivRnX* m_back = new_univ_rnx(params_glwe);
	PolyBiv* m_biv      = new_biv(params_glwe);
	rnx_random_vec(m, params_glwe);
	univ_rnx_to_biv(params_glwe, m_biv, m, 0);

	for (auto _ : state)
	{
		biv_to_univ_rnx(params_glwe, m_back, m_biv);
		benchmark::DoNotOptimize(m_back);
	}

	delete_univ_rnx(m);
	delete_univ_rnx(m_back);
	delete_biv(m_biv);

	delete_glwe_params(params_glwe);
	pvda_delete_backend(backend);
}

BENCHMARK(bench_biv_univ_rnx);
