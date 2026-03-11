#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "core/glwe/glwe_key.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  4
#define NLIMBSBASE (KBASE + 1) * 5
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -12

//! GLWE KEY PART (begin)

/**
 * @brief Ensures new_secret_key_values creates no NULL-pointer.
 */
Test(new_glwe_secret_key_values, basic)
{
	// Create a GLWE secret key's values
	PolyUniv** values = new_glwe_secret_key_values(NBASE, KBASE);

	// Assert values and its values is not NULL
	cr_assert(eq(int, values != NULL, 1));
	for (uint64_t j = 0; j < KBASE; j++) cr_assert(eq(int, values[j] != NULL, 1));

	// Clean up
	delete_glwe_secret_key_values(values, KBASE);
}

/**
 * @brief Ensures new_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
Test(new_glwe_secret_key, values_not_null)
{
	// Create a GLWE secret key
	GLWESecretKey* sk = new_glwe_secret_key(NBASE, KBASE);

	// Asserts the secret key is not NULL
	cr_assert(eq(int, sk != NULL, 1), "new_glwe_secret_key failed.");

	delete_glwe_secret_key(sk);
}

/**
 * @brief Ensures uniform_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
Test(uniform_glwe_secret_key, values_not_null)
{
	// Parameters
	MODULE* module = new_module_info(NBASE, FFT64);

	// Create a GLWE secret key
	GLWESecretKey* sk = new_glwe_secret_key(NBASE, KBASE);

	// Draw uniformly in Zn[X] the GLWE secret key's values
	int status = uniform_glwe_secret_key(module, sk, 2);

	// Asserts uniform_glwe_secret_key worked
	cr_assert(eq(int, status, 0), "uniform_glwe_secret_key failed.");

	// Clean up
	delete_module_info(module);
	delete_glwe_secret_key(sk);
}


//! GLWE KEY PART IN DFT SPACE (begin)

/**
 * @brief Ensure new_secret_key_values_dft creates no NULL-pointer.
 */
Test(new_glwe_secret_key_values_dft, basic)
{
	// Create a GLWE secret key's values
	PolyUnivDFT** values_dft = new_glwe_secret_key_values_dft(NBASE, KBASE);

	// Asserts values_dft is not NULL
	cr_assert(eq(int, values_dft != NULL, 1));
	for (uint64_t j = 0; j < KBASE; j++)
		cr_assert(eq(int, values_dft[j] != NULL, 1));

	// Clean up
	delete_glwe_secret_key_values_dft(values_dft, KBASE);
}

/**
 * @brief Ensures new_glwe_secret_key_dft returns a non-NULL pointer when values != NULL.
 */
Test(new_glwe_secret_key_dft, values_not_null)
{
	// Create a GLWE secret key in the DFT domain
	GLWESecretKeyDFT* sk_dft = new_glwe_secret_key_dft(NBASE, KBASE);

	// Asserts the secret key and its values are not NULL
	cr_assert(eq(int, sk_dft != NULL, 1));
	cr_assert(eq(int, sk_dft->values != NULL, 1));

	for (uint64_t j = 0; j < KBASE; j++) 
		cr_assert(eq(int, sk_dft->values[j] != NULL, 1));

	// Clean up
	delete_glwe_secret_key_dft(sk_dft);
}

/**
 * @brief Tests whether uniform_glwe_secret_key_gen works as intended.
 */
Test(uniform_glwe_secret_key_dft, what_s_inside)
{
	// Parameters
	MODULE* module = new_module_info(NBASE, FFT64);

	// Create a GLWE secret key 
	GLWESecretKeyDFT* sk_dft = new_glwe_secret_key_dft(NBASE, KBASE);

	// Draw uniformly the GLWE secret key values
	int status = uniform_glwe_secret_key_dft(module, sk_dft, 2);

	// Asserts uniform_glwe_secret_key_dft worked
	cr_assert(eq(int, status, 0), "uniform_glwe_secret_key_dft failed.");

	// Clean up
	delete_module_info(module);
	delete_glwe_secret_key_dft(sk_dft);
}
