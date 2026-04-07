#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdint.h>

#include "core/glwe/glwe_transform_key.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "univariate_polynomial.h"
#include "test_utils.h"

PvdaTstParams params = {1024, 1, 1, 1, 1, 0};

/**
 * @brief Tests whether transform_glwe_secret_key_not_dft_to_dft transforms the secret key out of the DFT domain, in DFT
 * space.
 */
Test(transform_glwe_secret_key_not_dft_to_dft, basic)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	GLWESecretKey* sk        = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft = alloc_glwe_secret_key_dft(params_glwe);
	PolyUniv* expected_poly  = new_univ(params_glwe);

	// Draws uniformly in Zn[X] the bivGLWE secret key's values
	uniform_glwe_secret_key(module, sk, 3);

	// Computes the bivGLWE secret key in the DFT domain
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	for (int k = 0; k < params_glwe->k; ++k)
	{
		univ_dft_to_coefs(module, expected_poly, glwe_sk_extract_poly_dft(sk_dft, k));
		for (int p = 0; p < params_glwe->nn; ++p)
		{
			cr_assert(eq(i64, expected_poly[p], glwe_sk_extract_poly(sk, k)[p]));
		}
	}

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);

	DELETE_PVDA_PARAMS_GLWE;
}
