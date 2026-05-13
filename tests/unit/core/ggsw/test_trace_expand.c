

#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <stdint.h>

#include "bivariate_polynomial.h"
#include "core/glwe/glwe_arithmetic.h"
#include "core/glwe/glwe_ciphertext.h"
#include "core/glwe/glwe_transform_key.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "glwegadget_arithmetic.h"
#include "glwegadget_key.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"
#include "utils.h"

PvdaParamTest(trace_expand, no_noise, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);

	PolyUnivRnX* m_univ_rnx     = new_univ_rnx(params_glwe);
	PolyUnivTnX* m_expected_tnx = new_univ_tnx(params_glwe);
	PolyUnivRnX* m_observed_rnx = new_univ_rnx(params_glwe);
	PolyBiv* biv_tmp            = new_biv_poly(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	int bundled[] = {2, 4, 3, 5, 14, params_glwe->nn, params_glwe->nn - 1};

	for (int i = 0; i < sizeof(bundled) / sizeof(bundled[0]); ++i)
	{
		int bund = bundled[i];

		memset(m_univ_rnx, 0, poly_univ_bytes(params_glwe));

		// Get the message in univariate RnX form for expected result
		normal_random_vec(m_univ_rnx, bund, 0, 0.01);

		glwe_secret_encrypt_rnx(module, glwe_ct, sk_prep, m_univ_rnx);

		GLWEAutomorphismKSK** ksks = calloc(2 * params_glwe->nn, sizeof(GLWEAutomorphismKSK*));
		for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
		{
			int64_t p = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
			ksks[p]   = new_automorphism_ksk(params_glwegadget);
			prepare_automorphism_key(module, ksks[p], sk_prep, p);
		}

		GLWECiphertext** results = calloc(bund, sizeof(GLWECiphertext*));
		for (int i = 0; i < bund; ++i)
		{
			results[i] = new_glwe(params_glwe);
		}

		glwegadget_trace_expand(module, results, bund, glwe_ct, (const GLWEAutomorphismKSK**)ksks, 2 * params_glwe->nn);

		int a = 1;
		for (int i = 0; i < bund; ++i)
		{
			glwe_secret_decrypt(module, biv_tmp, sk_prep, results[i]);
			biv_to_univ_rnx(params_glwe, m_observed_rnx, biv_tmp);
			int factor      = 1 << (32 - __builtin_clz(bund - 1));
			double expected = factor * m_univ_rnx[i];
			double actual   = m_observed_rnx[0];
			cr_assert(lt(dbl, rnx_torus_distance(expected, actual), 0.001));
			for (int p = 1; p < params_glwe->nn; ++p)
			{
				cr_assert(lt(dbl, rnx_torus_distance(0, m_observed_rnx[p]), 0.001));
			}
			delete_glwe(results[i]);
		}
		free(results);
		for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
		{
			int64_t p = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
			delete_automorphism_ksk(ksks[p]);
		}
		free(ksks);
	}

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	free(biv_tmp);
	delete_univ_rnx(m_univ_rnx);
	delete_univ_rnx(m_observed_rnx);
	delete_univ_tnx(m_expected_tnx);
	delete_glwe(glwe_ct);

	DELETE_PVDA_PARAMS_GGSWGAD;
}
