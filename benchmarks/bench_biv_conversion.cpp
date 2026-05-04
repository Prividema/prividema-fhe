#include <benchmark/benchmark.h>

#include <cmath>

extern "C" {
#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "rng.h"
#include "univariate_polynomial.h"
}

#define NBASE      (1 << 14)
#define KBASE      1
#define KAPPABASE  19
#define NLIMBSBASE (15 * 2)
#define LBASE      NLIMBSBASE / (KBASE + 1)

void test_univ_biv_rnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_FAST_UNIFORM);

	PolyUnivRnX* m = new_univ_rnx(params_glwe);
	PolyBiv* m_biv = new_biv_poly(params_glwe);
	normal_random_vec(m, NBASE, 0.0, 0.1);

	for (auto _ : state)
	{
		univ_rnx_to_biv(params_glwe, m_biv, m, 0);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_rnx(m);
	free(m_biv);

	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

BENCHMARK(test_univ_biv_rnx);

void test_biv_normalize(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_FAST_UNIFORM);

	PolyUnivRnX* m = new_univ_rnx(params_glwe);
	PolyBiv* m_biv = new_biv_poly(params_glwe);
	normal_random_vec(m, NBASE, 0.0, 0.1);

	univ_rnx_to_biv(params_glwe, m_biv, m, 0);

	for (auto _ : state)
	{
		pvda_vec_znx_normalize_base2k(module, params_glwe->kappa, m_biv, glwe_params_l_a(params_glwe), params_glwe->nn,
		                              m_biv, glwe_params_l_a(params_glwe), params_glwe->nn);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_rnx(m);
	free(m_biv);

	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
}

BENCHMARK(test_biv_normalize);

void test_univ_biv_rnx_via_tnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_FAST_UNIFORM);

	PolyUnivRnX* m     = new_univ_rnx(params_glwe);
	PolyBiv* m_biv     = new_biv_poly(params_glwe);
	PolyUnivTnX* m_tnx = new_univ_tnx(params_glwe);
	normal_random_vec(m, NBASE, 0.0, 0.1);

	for (auto _ : state)
	{
		univ_rnx_to_tnx(params_glwe, m_tnx, m);
		univ_tnx_to_biv(params_glwe, m_biv, m_tnx, 0);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_rnx(m);
	delete_univ_tnx(m_tnx);
	free(m_biv);

	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

BENCHMARK(test_univ_biv_rnx_via_tnx);

void test_univ_biv_tnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_FAST_UNIFORM);

	PolyUnivTnX* m = new_univ_tnx(params_glwe);
	PolyBiv* m_biv = new_biv_poly(params_glwe);
	uniform_random_pol_znx((PolyUniv*)m, NBASE, 64);

	for (auto _ : state)
	{
		univ_tnx_to_biv(params_glwe, m_biv, m, 0);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_tnx(m);
	free(m_biv);

	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

BENCHMARK(test_univ_biv_tnx);

void test_biv_univ_tnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_FAST_UNIFORM);

	PolyUnivTnX* m      = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_back = new_univ_tnx(params_glwe);
	PolyBiv* m_biv      = new_biv_poly(params_glwe);
	uniform_random_pol_znx((PolyUniv*)m, NBASE, 64);
	univ_tnx_to_biv(params_glwe, m_biv, m, 0);

	for (auto _ : state)
	{
		biv_to_univ_tnx(params_glwe, m_back, m_biv);
		benchmark::DoNotOptimize(m_back);
	}

	delete_univ_tnx(m_back);
	delete_univ_tnx(m);
	free(m_biv);

	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

BENCHMARK(test_biv_univ_tnx);

void test_biv_univ_rnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_FAST_UNIFORM);

	PolyUnivRnX* m      = new_univ_rnx(params_glwe);
	PolyUnivRnX* m_back = new_univ_rnx(params_glwe);
	PolyBiv* m_biv      = new_biv_poly(params_glwe);
	normal_random_vec(m, NBASE, 0.0, 0.1);
	univ_rnx_to_biv(params_glwe, m_biv, m, 0);

	for (auto _ : state)
	{
		biv_to_univ_rnx(params_glwe, m_back, m_biv);
		benchmark::DoNotOptimize(m_back);
	}

	delete_univ_rnx(m);
	delete_univ_rnx(m_back);
	free(m_biv);

	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

BENCHMARK(test_biv_univ_rnx);
