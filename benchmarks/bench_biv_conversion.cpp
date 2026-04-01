#include <benchmark/benchmark.h>

#include <cmath>

extern "C" {
#include "bivariate_polynomial.h"
#include "glwe.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
}

#define NBASE       (1 << 14)
#define KBASE       1
#define KAPPABASE   9
#define NLIMBSBASE  (15 * 2)
#define LBASE       NLIMBSBASE / (KBASE + 1)
#define SIGMABASE   -(LBASE / 2 + 1) * KAPPABASE

#define PROB_FACTOR 3

void test_univ_biv_rnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	PolyUnivRnX* m = new_univ_rnx(params_glwe);
	PolyBiv* m_biv = new_biv_poly(params_glwe);
	normal_random_vec(m, NBASE, 0.0, 0.1);

	for (auto _ : state)
	{
		univ_rnx_to_biv(params_glwe, m_biv, m);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_rnx(m);
	free(m_biv);
}

BENCHMARK(test_univ_biv_rnx);

void test_univ_biv_rnx_via_tnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	PolyUnivRnX* m     = new_univ_rnx(params_glwe);
	PolyBiv* m_biv     = new_biv_poly(params_glwe);
	PolyUnivTnX* m_tnx = new_univ_tnx(params_glwe);
	normal_random_vec(m, NBASE, 0.0, 0.1);

	for (auto _ : state)
	{
		univ_rnx_to_tnx(params_glwe, m_tnx, m);
		univ_tnx_to_biv(params_glwe, m_biv, m_tnx);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_rnx(m);
	free(m_biv);
}

BENCHMARK(test_univ_biv_rnx_via_tnx);

void test_univ_biv_tnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	PolyUnivTnX* m = new_univ_tnx(params_glwe);
	PolyBiv* m_biv = new_biv_poly(params_glwe);
	uniform_random_pol_znx((PolyUniv*)m, NBASE, 64);

	for (auto _ : state)
	{
		univ_tnx_to_biv(params_glwe, m_biv, m);
		benchmark::DoNotOptimize(m_biv);
	}

	delete_univ_tnx(m);
	free(m_biv);
}

BENCHMARK(test_univ_biv_tnx);

void test_biv_univ_tnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	PolyUnivTnX* m      = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_back = new_univ_tnx(params_glwe);
	PolyBiv* m_biv      = new_biv_poly(params_glwe);
	uniform_random_pol_znx((PolyUniv*)m, NBASE, 64);
	univ_tnx_to_biv(params_glwe, m_biv, m);

	for (auto _ : state)
	{
		biv_to_univ_tnx(params_glwe, m_back, m_biv);
		benchmark::DoNotOptimize(m_back);
	}

	delete_univ_tnx(m_back);
	delete_univ_tnx(m);
	free(m_biv);
}

BENCHMARK(test_biv_univ_tnx);

void test_biv_univ_rnx(benchmark::State& state)
{
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	PolyUnivRnX* m      = new_univ_rnx(params_glwe);
	PolyUnivRnX* m_back = new_univ_rnx(params_glwe);
	PolyBiv* m_biv      = new_biv_poly(params_glwe);
	normal_random_vec(m, NBASE, 0.0, 0.1);
	univ_rnx_to_biv(params_glwe, m_biv, m);

	for (auto _ : state)
	{
		biv_to_univ_rnx(params_glwe, m_back, m_biv);
		benchmark::DoNotOptimize(m_back);
	}

	delete_univ_rnx(m);
	delete_univ_rnx(m_back);
	free(m_biv);
}

BENCHMARK(test_biv_univ_rnx);
