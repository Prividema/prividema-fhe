#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "core/glwe/glwe_key.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  4
#define NLIMBSBASE (KBASE + 1) * 5
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -12

//! bivGLWE KEY PART (begin)

/**
 * @brief Ensures new_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
Test(new_glwe_secret_key, values_not_null)
{
	GLWESecretKey* sk = alloc_glwe_secret_key(NBASE, KBASE);

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
	MODULE* module = pvda_new_module_info(NBASE);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = alloc_glwe_secret_key(NBASE, KBASE);

	// Draw uniformly in Zn[X] the bivGLWE secret key's values
	int status = uniform_glwe_secret_key(module, sk, 2);

	// Asserts uniform_glwe_secret_key worked
	cr_assert(eq(int, status, 0), "uniform_glwe_secret_key failed.");

	// Clean up
	pvda_delete_module_info(module);
	delete_glwe_secret_key(sk);
}
