#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "core/glwe/glwe_transform_key.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  4
#define NLIMBSBASE (KBASE + 1) * 5
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -12

//! GLWE KEY PART (begin)

/**
 * @brief Tests wether transform_glwe_secret_key_dft_to_not_dft transforms the secret key in the DFT domain, out of DFT
 * space.
 */
Test(transform_glwe_secret_key_dft_to_not_dft, basic)
{
	// Parameters
	MODULE* module = new_module_info(NBASE, FFT64);

	// Create a GLWE secret key 
	GLWESecretKeyDFT* sk_dft = new_glwe_secret_key_dft(NBASE, KBASE);
	GLWESecretKey* sk =	new_glwe_secret_key(NBASE, KBASE);

	// Draws uniformly in Zn[X] the secret key's values
	uniform_glwe_secret_key_dft(module, sk_dft, 3);

	// Compute the GLWE secret key out of DFT space
	int status = transform_glwe_secret_key_dft_to_not_dft(module, sk, sk_dft);

	// Asserts transform_glwe_secret_key_dft_to_not_dft succeed
	cr_assert(eq(int, status, 0), "transform_glwe_secret_key_dft_to_not_dft failed.");

	// Clean up
	delete_module_info(module);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_secret_key(sk);
}


//! GLWE KEY PART IN DFT SPACE (begin)

/**
 * @brief Tests wether transform_glwe_secret_key_not_dft_to_dft transforms the secret key out of the DFT domain, in DFT
 * space.
 */
Test(transform_glwe_secret_key_not_dft_to_dft, basic)
{
	// Parameters
	MODULE* module = new_module_info(NBASE, FFT64);

	// Variables
	GLWESecretKey* sk = new_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft = new_glwe_secret_key_dft(NBASE, KBASE);
	
	// Draws uniformly in Zn[X] the GLWE secret key's values
	uniform_glwe_secret_key(module, sk, 3);

	// Computes the GLWE secret key in the DFT domain
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// Clean up
	delete_module_info(module);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
}
