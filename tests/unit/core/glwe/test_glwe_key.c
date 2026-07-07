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

/**
 * @brief Ensures uniform_glwe_secret_key's coefficients are in the provided bounds
 */
PvdaParamTest(uniform_glwe_secret_key, values_bounds, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = alloc_glwe_secret_key(params_glwe);

	int bitsp[] = {1, 2, 3, 4, 5, 6, 7, 11, 13, 32, 33};

	for (int i = 0; i < sizeof(bitsp) / sizeof(bitsp[0]); ++i)
	{
		int nb_bits = bitsp[i];
		int status  = uniform_glwe_secret_key(module, sk, nb_bits);

		int64_t lb = -(1 << (nb_bits - 1));
		int64_t ub = -lb - 1;

		cr_assert(eq(int, status, 0), "uniform_glwe_secret_key failed.");

		for (int i = 0; i < params_glwe->nn * params_glwe->k; ++i)
		{
			cr_assert(sk->values[i] >= lb || sk->values[i] <= ub, "out-of-bounds values");
		}
	}

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Ensures binary_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
PvdaParamTest(binary_glwe_secret_key, values_not_null, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = alloc_glwe_secret_key(params_glwe);

	// Draw uniformly in Zn[X] the bivGLWE secret key's values
	int status = binary_glwe_secret_key(module, sk);

	// Asserts uniform_glwe_secret_key worked
	cr_assert(eq(int, status, 0), "uniform_glwe_secret_key failed.");

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Ensures binary_glwe_secret_key's generated key is actually binary
 */
PvdaParamTest(binary_glwe_secret_key, values_binary, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = alloc_glwe_secret_key(params_glwe);

	int status = binary_glwe_secret_key(module, sk);

	cr_assert(eq(int, status, 0), "uniform_glwe_secret_key failed.");

	for (int i = 0; i < params_glwe->nn * params_glwe->k; ++i)
	{
		cr_assert(sk->values[i] == 0 || sk->values[i] == 1, "Non-binary values");
	}

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}
