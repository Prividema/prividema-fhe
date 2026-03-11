#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "core/ggsw/ggsw_transform_key.h"

#define NBASE            1024
#define KBASE            1
#define KAPPABASE        4
#define NLIMBSBASE       (KBASE + 1) * 2
#define LBASE            NLIMBSBASE / (KBASE + 1)
#define SIGMABASE        -(LBASE / 2 + 1) * KAPPABASE

#define K_TILDEBASE      1
#define KAPPA_TILDEBASE  4
#define NLIMBS_TILDEBASE (K_TILDEBASE + 1) * 2
#define L_TILDEBASE      NLIMBS_TILDEBASE / (K_TILDEBASE + 1)
#define SIGMA_TILDEBASE  -3


//! GGSW KEY PART (begin)


/**
 * @brief Tests wether transform_ggsw_secret_key_dft_to_not_dft transforms the secret key in the DFT domain, out of DFT
 * space.
 */
Test(transform_ggsw_secret_key_dft_to_not_dft, basic)
{
	// Parameters
	MODULE* module = new_module_info(NBASE, FFT64);

	// Create a GGSW secret key 
	GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NBASE, KBASE);
	GGSWSecretKey* sk =	new_ggsw_secret_key(NBASE, KBASE);

	// Draw uniformly in Zn[X] the secret key's values
	uniform_ggsw_secret_key_dft(module, sk_dft, 3);

	// Compute the GGSW secret key out of DFT space
	int status = transform_ggsw_secret_key_dft_to_not_dft(module, sk, sk_dft);

	// Asserts transform_ggsw_secret_key_dft_to_not_dft succeed
	cr_assert(eq(int, status, 0), "transform_ggsw_secret_key_dft_to_not_dft failed.");

	// Clean up
	delete_module_info(module);
	delete_ggsw_secret_key_dft(sk_dft);
	delete_ggsw_secret_key(sk);
}

/**
 * @brief Tests wether transform_ggsw_secret_key_to_glwe_secret_key actually runs.
 * 
 */
Test(transform_ggsw_secret_key_to_glwe_secret_key, basic)
{
	//Parameters
	MODULE* module = new_module_info(NBASE, FFT64);

	// Variables
	GGSWSecretKey* sk_ggsw = new_ggsw_secret_key(NBASE, KBASE);
	GLWESecretKey* sk_glwe = new_glwe_secret_key(NBASE, KBASE);

	// Draws uniformly in Zn[X] the GGSW secret key values
	uniform_ggsw_secret_key(module, sk_ggsw, 2);

	// Computes the result GLWE secret key
	transform_ggsw_secret_key_to_glwe_secret_key(sk_glwe, sk_ggsw);

	// Clean up
	delete_glwe_secret_key(sk_glwe);
	delete_ggsw_secret_key(sk_ggsw);
	delete_module_info(module);
}


//! GGSW KEY PART IN DFT SPACE (begin)

/**
 * @brief Tests wether transform_ggsw_secret_key_not_dft_to_dft transforms the secret key out of the DFT domain, in DFT
 * space.
 */
Test(transform_ggsw_secret_key_not_dft_to_dft, basic)
{
	// Parameters
	MODULE* module = new_module_info(NBASE, FFT64);

	// Variables
	GGSWSecretKey* sk = new_ggsw_secret_key(NBASE, KBASE);
	GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NBASE, KBASE);
	
	// Draws uniformly in Zn[X] the GGSW secret key's values
	uniform_ggsw_secret_key(module, sk, 3);

	// Computes the GGSW secret key in the DFT domain
	transform_ggsw_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// Clean up
	delete_module_info(module);
	delete_ggsw_secret_key(sk);
	delete_ggsw_secret_key_dft(sk_dft);
}

/**
 * @brief Tests wether transform_ggsw_secret_key_dft_to_glwe_secret_key_dft actually runs.
 * 
 */
Test(transform_ggsw_secret_key_dft_to_glwe_secret_key_dft, basic)
{
	//Parameters
	MODULE* module = new_module_info(NBASE, FFT64);

	// Variables
	GGSWSecretKeyDFT* sk_ggsw_dft = new_ggsw_secret_key_dft(NBASE, KBASE);
	GLWESecretKeyDFT* sk_glwe_dft = new_glwe_secret_key_dft(NBASE, KBASE);

	// Draws uniformly in Zn[X] the GGSW secret key values in the DFT domain
	uniform_ggsw_secret_key_dft(module, sk_ggsw_dft, 2);

	// Computes the result GLWE secret key in the DFT domain
	transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(sk_glwe_dft, sk_ggsw_dft);

	// Clean up
	delete_glwe_secret_key_dft(sk_glwe_dft);
	delete_ggsw_secret_key_dft(sk_ggsw_dft);
	delete_module_info(module);
}