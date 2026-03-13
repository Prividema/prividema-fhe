#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "core/ggsw/ggsw_key.h"

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

//! bivGGSW KEY PART (begin)

/**
 * @brief Ensures new_secret_key_values creates no NULL-pointer.
 */
Test(new_ggsw_secret_key_values, basic)
{
	// Create a bivGGSW secret key's values
	PolyUniv** values = new_ggsw_secret_key_values(NBASE, KBASE);

	// Assert values and its values is not NULL
	cr_assert(eq(int, values != NULL, 1));
	for (uint64_t j = 0; j < KBASE; j++) cr_assert(eq(int, values[j] != NULL, 1));

	// Clean up
	delete_ggsw_secret_key_values(values, KBASE);
}

/**
 * @brief Ensures new_ggsw_secret_key returns a non-NULL pointer when values != NULL.
 */
Test(new_ggsw_secret_key, values_not_null)
{
	// Create a bivGGSW secret key
	GGSWSecretKey* sk = new_ggsw_secret_key(NBASE, KBASE);

	// Asserts the secret key is not NULL
	cr_assert(eq(int, sk != NULL, 1), "new_ggsw_secret_key failed.");

	delete_ggsw_secret_key(sk);
}

/**
 * @brief Ensures uniform_ggsw_secret_key returns a non-NULL pointer when values != NULL.
 */
Test(uniform_ggsw_secret_key, values_not_null)
{
	// Parameters
	MODULE* module = new_module_info(NBASE, FFT64);

	// Create a bivGGSW secret key
	GGSWSecretKey* sk = new_ggsw_secret_key(NBASE, KBASE);

	// Draw uniformly in Zn[X] the bivGGSW secret key's values
	int status = uniform_ggsw_secret_key(module, sk, 2);

	// Asserts uniform_ggsw_secret_key worked
	cr_assert(eq(int, status, 0), "uniform_ggsw_secret_key failed.");

	// Clean up
	delete_module_info(module);
	delete_ggsw_secret_key(sk);
}


//! bivGGSW KEY PART IN DFT SPACE (begin)

/**
 * @brief Ensure new_secret_key_values_dft creates no NULL-pointer.
 */
Test(new_ggsw_secret_key_values_dft, basic)
{
	// Create a bivGGSW secret key's values
	PolyUnivDFT** values_dft = new_ggsw_secret_key_values_dft(NBASE, KBASE);

	// Asserts values_dft is not NULL
	cr_assert(eq(int, values_dft != NULL, 1));
	for (uint64_t j = 0; j < KBASE; j++)
		cr_assert(eq(int, values_dft[j] != NULL, 1));

	// Clean up
	delete_ggsw_secret_key_values_dft(values_dft, KBASE);
}

/**
 * @brief Ensures new_ggsw_secret_key_dft returns a non-NULL pointer when values != NULL.
 */
Test(new_ggsw_secret_key_dft, values_not_null)
{
	// Create a bivGGSW secret key in the DFT domain
	GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NBASE, KBASE);

	// Asserts the secret key and its values are not NULL
	cr_assert(eq(int, sk_dft != NULL, 1));
	cr_assert(eq(int, sk_dft->values != NULL, 1));

	for (uint64_t j = 0; j < KBASE; j++) 
		cr_assert(eq(int, sk_dft->values[j] != NULL, 1));

	// Clean up
	delete_ggsw_secret_key_dft(sk_dft);
}

/**
 * @brief Tests whether uniform_ggsw_secret_key_gen works as intended.
 */
Test(uniform_ggsw_secret_key_dft, what_s_inside)
{
	// Parameters
	MODULE* module = new_module_info(NBASE, FFT64);

	// Create a bivGGSW secret key 
	GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NBASE, KBASE);

	// Draw uniformly the bivGGSW secret key values
	int status = uniform_ggsw_secret_key_dft(module, sk_dft, 2);

	// Asserts uniform_ggsw_secret_key_dft worked
	cr_assert(eq(int, status, 0), "uniform_ggsw_secret_key_dft failed.");

	// Clean up
	delete_module_info(module);
	delete_ggsw_secret_key_dft(sk_dft);
}
