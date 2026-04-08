#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "core/glwe/glwe_key.h"
#include "test_utils.h"

//! bivGLWE KEY PART (begin)

/**
 * @brief Ensures new_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
PvdaParamTest(new_glwe_secret_key, values_not_null, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	GLWESecretKey* sk = alloc_glwe_secret_key(params_glwe);

	cr_assert(eq(int, sk != NULL, 1), "new_glwe_secret_key failed.");

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Ensures uniform_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
PvdaParamTest(uniform_glwe_secret_key, values_not_null, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = alloc_glwe_secret_key(params_glwe);

	// Draw uniformly in Zn[X] the bivGLWE secret key's values
	int status = uniform_glwe_secret_key(module, sk, 2);

	// Asserts uniform_glwe_secret_key worked
	cr_assert(eq(int, status, 0), "uniform_glwe_secret_key failed.");

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}
