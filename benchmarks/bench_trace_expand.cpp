
#include <benchmark/benchmark.h>

#include <cmath>
#include <cstring>

#include "utils.hpp"

extern "C" {
#include "backend.h"
#include "bivariate_polynomial.h"
#include "ggsw_arithmetic.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "glwegadget_arithmetic.h"
#include "glwegadget_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"
}

#include "params.h"

#define D 16  //Number of non-zero coefficients in the GLWE to expand

void test_expand_trace(benchmark::State& state)
{
	PvdaBackend* module = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* params_glwegadget = new_glwegadget_params(params_glwe, KAPPABASE, LBASE);
	params_glwe->fast_uniform_nb_bits   = 0;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);

	PolyUnivRnX* m_univ_rnx = new_univ_rnx(params_glwe);
	GLWECiphertext* glwe_ct = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	GLWEAutomorphismKeyCollection* auto_key_collection = new_automorphism_key_collection(2 * params_glwe->nn);
	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                     = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* auto_key = new_automorphism_key(params_glwegadget);
		compute_automorphism_key(module, auto_key, sk_prep, p);
		glwegadget_key_collection_put_key(auto_key_collection, auto_key, p);
	}
	int bund = D;

	memset(m_univ_rnx, 0, poly_univ_bytes(params_glwe));

	rnx_random_vec(module, m_univ_rnx, params_glwe);

	glwe_secret_encrypt_rnx(module, glwe_ct, sk_prep, m_univ_rnx);

	GLWECiphertext** results = (GLWECiphertext**)calloc(bund, sizeof(GLWECiphertext*));
	for (int i = 0; i < bund; ++i)
	{
		results[i] = new_glwe(params_glwe);
	}

	for (auto _ : state)
	{
		glwe_trace_expand(module, results, bund, glwe_ct, auto_key_collection);

		benchmark::DoNotOptimize(results);
	}

	for (int i = 0; i < bund; ++i)
	{
		delete_glwe(results[i]);
	}
	free(results);

	delete_automorphism_key_collection(auto_key_collection, 1);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	delete_univ_rnx(m_univ_rnx);
	delete_glwe(glwe_ct);

	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_glwegadget_params(params_glwegadget);
}

BENCHMARK(test_expand_trace);

void test_hom_trace(benchmark::State& state)
{
	PvdaBackend* module = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* params_glwegadget = new_glwegadget_params(params_glwe, KAPPABASE, LBASE);
	params_glwe->fast_uniform_nb_bits   = 0;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);

	PolyUnivRnX* m_univ_rnx = new_univ_rnx(params_glwe);
	GLWECiphertext* glwe_ct = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	GLWEAutomorphismKeyCollection* auto_key_collection = new_automorphism_key_collection(2 * params_glwe->nn);
	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                     = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* auto_key = new_automorphism_key(params_glwegadget);
		compute_automorphism_key(module, auto_key, sk_prep, p);
		glwegadget_key_collection_put_key(auto_key_collection, auto_key, p);
	}
	int bund = D;

	memset(m_univ_rnx, 0, poly_univ_bytes(params_glwe));

	rnx_random_vec(module, m_univ_rnx, params_glwe);

	glwe_secret_encrypt_rnx(module, glwe_ct, sk_prep, m_univ_rnx);

	GLWECiphertext* result;
	result = new_glwe(params_glwe);

	for (auto _ : state)
	{
		glwe_hom_trace(module, result, bund, D / 2, glwe_ct, auto_key_collection);

		benchmark::DoNotOptimize(result);
	}

	delete_glwe(result);

	delete_automorphism_key_collection(auto_key_collection, 1);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	delete_univ_rnx(m_univ_rnx);
	delete_glwe(glwe_ct);

	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_glwegadget_params(params_glwegadget);
}

BENCHMARK(test_hom_trace);

void test_expand_compressed_trace(benchmark::State& state)
{
	PvdaBackend* module = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* params_glwegadget = new_glwegadget_params(params_glwe, KAPPABASE, LBASE);
	GGSWParams* params_ggsw             = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NLIMBSBASE);

	uint64_t k  = params_glwe->k;
	uint64_t nn = params_glwe->nn;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);

	PolyUniv* m_univ        = new_univ(params_glwe);
	PolyBiv* biv_tmp        = new_biv(params_glwe);
	GLWECiphertext* glwe_ct = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	GLWEAutomorphismKeyCollection* automorphism_keys = new_automorphism_key_collection(2 * params_glwe->nn);

	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                     = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* auto_key = new_automorphism_key(params_glwegadget);
		compute_automorphism_key(module, auto_key, sk_prep, p);
		glwegadget_key_collection_put_key(automorphism_keys, auto_key, p);
	}

	GGSWCiphertextPrep** ggsw_ksks = (GGSWCiphertextPrep**)calloc(k, sizeof(GGSWCiphertextPrep*));
	GGSWCiphertext* ggsw_tmp       = new_ggsw(params_ggsw);
	PolyUniv* neg_sk_i             = new_univ(params_glwe);
	for (uint64_t i = 0; i < k; ++i)
	{
		ggsw_ksks[i] = new_ggsw_prep(params_ggsw);
		for (int p = 0; p < nn; ++p)
		{
			neg_sk_i[p] = -glwe_prepared_sk_extract_poly_coefs(sk_prep, i)[p];
		}
		ggsw_secret_encrypt(module, ggsw_tmp, sk_prep, neg_sk_i);
		ggsw_prepare(module, ggsw_ksks[i], ggsw_tmp);
	}
	delete_univ(neg_sk_i);
	delete_ggsw(ggsw_tmp);

	int bund = D;
	memset(m_univ, 0, poly_univ_bytes(params_glwe));

	uniform_pow2_random_pol_znx(module, m_univ, bund, 1);

	glwegadget_packed_secret_encrypt(module, glwe_ct, params_glwegadget, sk_prep, m_univ, bund);

	GGSWCiphertext** results = (GGSWCiphertext**)calloc(bund, sizeof(GGSWCiphertext*));
	for (int i = 0; i < bund; ++i)
	{
		results[i] = new_ggsw(params_ggsw);
	}

	for (auto _ : state)
	{
		packed_glwegadget_trace_expand_ggsw(module, results, bund, params_glwegadget->l_tilde, glwe_ct,
		                                    automorphism_keys, (const GGSWCiphertextPrep**)ggsw_ksks);
		benchmark::DoNotOptimize(results);
	}
	for (int i = 0; i < bund; ++i)
	{
		delete_ggsw(results[i]);
	}
	delete_automorphism_key_collection(automorphism_keys, 1);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	delete_biv(biv_tmp);
	delete_glwe(glwe_ct);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_glwegadget_params(params_glwegadget);
	delete_ggsw_params(params_ggsw);
}
BENCHMARK(test_expand_compressed_trace);

void test_expand_compressed_trace_gad(benchmark::State& state)
{
	PvdaBackend* module = pvda_new_spqlios_backend(NBASE);
	GLWEParams* params_glwe =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* params_glwegadget = new_glwegadget_params(params_glwe, KAPPABASE, LBASE);
	GGSWParams* params_ggsw             = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NLIMBSBASE);

	uint64_t k  = params_glwe->k;
	uint64_t nn = params_glwe->nn;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);

	PolyUniv* m_univ        = new_univ(params_glwe);
	PolyBiv* biv_tmp        = new_biv(params_glwe);
	GLWECiphertext* glwe_ct = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	GLWEAutomorphismKeyCollection* auto_key_collection = new_automorphism_key_collection(2 * params_glwe->nn);

	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                     = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* auto_key = new_automorphism_key(params_glwegadget);
		compute_automorphism_key(module, auto_key, sk_prep, p);
		glwegadget_key_collection_put_key(auto_key_collection, auto_key, p);
	}

	GGSWCiphertext** ggsw_ksks = (GGSWCiphertext**)calloc(k, sizeof(GGSWCiphertext*));
	PolyUniv* neg_sk_i         = new_univ(params_glwe);
	for (uint64_t i = 0; i < k; ++i)
	{
		ggsw_ksks[i] = new_ggsw(params_ggsw);
		for (int p = 0; p < nn; ++p)
		{
			neg_sk_i[p] = -glwe_prepared_sk_extract_poly_coefs(sk_prep, i)[p];
		}
		ggsw_secret_encrypt(module, ggsw_ksks[i], sk_prep, neg_sk_i);
	}
	delete_univ(neg_sk_i);

	int bund = D;

	memset(m_univ, 0, poly_univ_bytes(params_glwe));

	// Get the message in univariate RnX form for expected result
	uniform_pow2_random_pol_znx(module, m_univ, bund, 1);

	glwegadget_packed_secret_encrypt(module, glwe_ct, params_glwegadget, sk_prep, m_univ, bund);

	GLWEGadgetCiphertext** results = (GLWEGadgetCiphertext**)calloc(bund, sizeof(GGSWCiphertext*));
	for (int i = 0; i < bund; ++i)
	{
		results[i] = new_glwegadget(params_glwegadget);
	}

	for (auto _ : state)
	{
		packed_glwegadget_trace_expand(module, results, bund, params_glwegadget->l_tilde, glwe_ct, auto_key_collection);

		benchmark::DoNotOptimize(results);
	}
	for (int i = 0; i < bund; ++i)
	{
		delete_glwegadget(results[i]);
	}
	delete_automorphism_key_collection(auto_key_collection, 1);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	delete_biv(biv_tmp);
	delete_glwe(glwe_ct);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_glwegadget_params(params_glwegadget);
	delete_ggsw_params(params_ggsw);
}
BENCHMARK(test_expand_compressed_trace_gad);
