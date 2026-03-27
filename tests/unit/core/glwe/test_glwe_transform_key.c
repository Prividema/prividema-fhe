#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdint.h>

#include "core/glwe/glwe_transform_key.h"
#include "glwe_key.h"
#include "univariate_polynomial.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  4
#define NLIMBSBASE (KBASE + 1) * 5
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -12

/**
 * @brief Tests wether transform_glwe_secret_key_not_dft_to_dft transforms the secret key out of the DFT domain, in DFT
 * space.
 */
Test(transform_glwe_secret_key_not_dft_to_dft, basic)
{
	// Parameters
	MODULE* module = pvda_new_module_info(NBASE);

	// Variables
	GLWESecretKey* sk        = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft = alloc_glwe_secret_key_dft(NBASE, KBASE);
	PolyUniv* expected_poly  = calloc(NBASE, sizeof(int64_t));  // TODO: correct type

	// Draws uniformly in Zn[X] the bivGLWE secret key's values
	uniform_glwe_secret_key(module, sk, 3);

	// Computes the bivGLWE secret key in the DFT domain
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	for (int k = 0; k < KBASE; ++k)
	{
		univ_dft_to_coefs(module, expected_poly, glwe_sk_extract_poly_dft(sk_dft, k));
		for (int p = 0; p < NBASE; ++p)
		{
			cr_assert(eq(i64, expected_poly[p], glwe_sk_extract_poly(sk, k)[p]));
		}
	}

	// Clean up
	pvda_delete_module_info(module);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
}
