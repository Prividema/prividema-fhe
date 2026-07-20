#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdint.h>

#include "core/glwe/glwe_transform_key.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

PvdaTstParams params = {1024, 1, 1, 1, 1, 0};

/**
 * @brief Tests whether transform_glwe_secret_key_not_dft_to_dft transforms the secret key out of the DFT domain, in DFT
 * space.
 */
PvdaParamTest(transform_glwe_secret_key_not_dft_to_dft, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	PolyUniv* expected_poly        = new_univ(params_glwe);

	// Draws uniformly in Zn[X] the bivGLWE secret key's values
	uniform_glwe_secret_key(module, sk, 3);

	// Computes the bivGLWE secret key in the DFT domain
	glwe_sk_prepare(module, sk_prep, sk);

	for (int k = 0; k < params_glwe->k; ++k)
	{
		univ_dft_to_coefs(module, expected_poly, glwe_prepared_sk_extract_poly_dft(sk_prep, k));
		for (int p = 0; p < params_glwe->nn; ++p)
		{
			cr_assert(eq(i64, expected_poly[p], glwe_sk_extract_poly(sk, k)[p]));
		}
		for (int p = 0; p < params_glwe->nn; ++p)
		{
			cr_assert(eq(i64, glwe_prepared_sk_extract_poly_coefs(sk_prep, k)[p], glwe_sk_extract_poly(sk, k)[p]));
		}
	}

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_univ(expected_poly);

	DELETE_PVDA_PARAMS_GLWE;
}
