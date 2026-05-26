#include <benchmark/benchmark.h>

#include <cmath>

extern "C" {
#include "bivariate_polynomial.h"
#include "ggsw_arithmetic.h"
#include "ggsw_params.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "glwegadget_arithmetic.h"
#include "glwegadget_ciphertext.h"
#include "rng.h"
#include "univariate_polynomial.h"
}

#include "params.h"

void test_glwegad_half_prod(benchmark::State& state)
{
	MODULE* module = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);
	GGSWParams* params_ggsw             = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NLIMBSBASE);
	GLWEGadgetParams* params_glwegadget = new_glwegadget_params(params_glwe, KAPPABASE, LBASE);

	GLWESecretKey* sk                      = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft               = alloc_glwe_secret_key_dft(params_glwe);
	GLWECiphertext* glwe                   = new_glwe(params_glwe);
	GLWEGadgetCiphertext* glwegad          = new_glwegadget(params_glwegadget);
	GLWEGadgetCiphertextPrep* glwegad_prep = new_glwegadget_prep(params_glwegadget);
	PolyUniv* u_univ                       = new_univ(params_glwe);
	PolyUnivTnX* m_univ_tnx                = new_univ_tnx(params_glwe);
	PolyBiv* m                             = new_biv_poly(params_glwe);

	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);
	uniform_random_pol_znx(u_univ, params_glwe->nn, 3);
	uniform_random_pol_znx((PolyUniv*)m_univ_tnx, params_glwe->nn, 12);
	univ_tnx_to_biv(params_glwe, m, m_univ_tnx, 0);

	glwegadget_secret_encrypt(module, glwegad, sk_dft, u_univ);
	glwegadget_prepare(module, glwegad_prep, glwegad);

	for (auto _ : state)
	{
		glwegadget_half_prod(module, glwe, glwegad_prep, m);
		benchmark::DoNotOptimize(glwe);
	}

	free(m);
	delete_univ(u_univ);
	delete_univ_tnx(m_univ_tnx);

	delete_glwegadget(glwegad);
	delete_glwegadget_prep(glwegad_prep);
	delete_glwe(glwe);

	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
	delete_glwegadget_params(params_glwegadget);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
}

BENCHMARK(test_glwegad_half_prod);

void test_glwegad_half_prod_dft(benchmark::State& state)
{
	double sigma = ldexp(1.0, 4 - (LBASE)*KAPPABASE);

	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NLIMBSBASE);
	GLWEGadgetParams* params_glwegadget = new_glwegadget_params(params_glwe, KAPPABASE, LBASE);

	GLWESecretKey* sk                      = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft               = alloc_glwe_secret_key_dft(params_glwe);
	GLWECiphertext* glwe                   = new_glwe(params_glwe);
	GLWECiphertextDFT* glwe_dft            = new_glwe_dft(params_glwe);
	GLWEGadgetCiphertext* glwegad          = new_glwegadget(params_glwegadget);
	GLWEGadgetCiphertextPrep* glwegad_prep = new_glwegadget_prep(params_glwegadget);
	PolyUniv* u_univ                       = new_univ(params_glwe);
	PolyUnivTnX* m_univ_tnx                = new_univ_tnx(params_glwe);
	PolyBiv* m                             = new_biv_poly(params_glwe);
	PolyBivDFT* m_dft                      = new_biv_poly_dft(params_glwe);

	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);
	uniform_random_pol_znx(u_univ, params_glwe->nn, 3);
	uniform_random_pol_znx((PolyUniv*)m_univ_tnx, params_glwe->nn, 12);
	univ_tnx_to_biv(params_glwe, m, m_univ_tnx, 0);
	biv_coefs_to_dft(module, params_glwe, m_dft, m);

	glwegadget_secret_encrypt(module, glwegad, sk_dft, u_univ);
	glwegadget_prepare(module, glwegad_prep, glwegad);

	for (auto _ : state)
	{
		glwegadget_half_prod_dft_to_dft(module, glwe_dft, glwegad_prep, m_dft);
		benchmark::DoNotOptimize(glwe);
	}

	free(m);
	free(m_dft);
	delete_univ(u_univ);
	delete_univ_tnx(m_univ_tnx);

	delete_glwegadget(glwegad);
	delete_glwegadget_prep(glwegad_prep);
	delete_glwe(glwe);
	delete_glwe_dft(glwe_dft);

	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
	delete_glwegadget_params(params_glwegadget);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
}

BENCHMARK(test_glwegad_half_prod_dft);
