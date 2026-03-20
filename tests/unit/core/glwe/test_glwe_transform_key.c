#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "core/glwe/glwe_transform_key.h"

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
	MODULE* module = new_module_info(NBASE, FFT64);

	// Variables
	GLWESecretKey* sk        = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft = alloc_glwe_secret_key_dft(NBASE, KBASE);

	// Draws uniformly in Zn[X] the bivGLWE secret key's values
	uniform_glwe_secret_key(module, sk, 3);

	// Computes the bivGLWE secret key in the DFT domain
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// Clean up
	delete_module_info(module);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
}
