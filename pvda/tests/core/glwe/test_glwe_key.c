#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "backend_arithmetic.h"
#include "bivariate_polynomial.h"
#include "core/glwe/glwe_key.h"
#include "logger.h"
#include "test_utils.h"

//! bivGLWE KEY PART (begin)

/**
 * @brief Ensures new_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
PvdaParamTest(new_glwe_secret_key, values_not_null, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	GLWESecretKey* sk = new_glwe_secret_key(params_glwe);

	cr_assert(eq(int, sk != NULL, 1), "new_glwe_secret_key failed.");

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Ensures glwe_sk_random_uniform returns a non-NULL pointer when values != NULL.
 */
PvdaParamTest(glwe_sk_random_uniform, values_not_null, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = new_glwe_secret_key(params_glwe);

	// Draw uniformly in Zn[X] the bivGLWE secret key's values
	int status = glwe_sk_random_uniform(module, sk, 2);

	// Asserts glwe_sk_random_uniform worked
	cr_assert(eq(int, status, 0), "glwe_sk_random_uniform failed.");

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Ensures glwe_sk_random_uniform's coefficients are in the provided bounds
 */
PvdaParamTest(glwe_sk_random_uniform, values_bounds, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = new_glwe_secret_key(params_glwe);

	int bitsp[] = {1, 2, 3, 4, 5, 6, 7, 11, 13, 32, 33};

	for (int i = 0; i < sizeof(bitsp) / sizeof(bitsp[0]); ++i)
	{
		int nb_bits = bitsp[i];
		int status  = glwe_sk_random_uniform(module, sk, nb_bits);

		int64_t lb = -(1 << (nb_bits - 1));
		int64_t ub = -lb - 1;

		cr_assert(eq(int, status, 0), "glwe_sk_random_uniform failed.");

		for (int i = 0; i < params_glwe->nn * params_glwe->k; ++i)
		{
			cr_assert(sk->values[i] >= lb || sk->values[i] <= ub, "out-of-bounds values");
		}
	}

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Ensures glwe_sk_random_binary returns a non-NULL pointer when values != NULL.
 */
PvdaParamTest(glwe_sk_random_binary, values_not_null, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = new_glwe_secret_key(params_glwe);

	// Draw uniformly in Zn[X] the bivGLWE secret key's values
	int status = glwe_sk_random_binary(module, sk);

	cr_assert(eq(int, status, 0), "glwe_sk_random_binary failed.");

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Ensures glwe_sk_random_binary's generated key is actually binary
 */
PvdaParamTest(glwe_sk_random_binary, values_binary, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = new_glwe_secret_key(params_glwe);

	int status = glwe_sk_random_binary(module, sk);

	cr_assert(eq(int, status, 0), "glwe_sk_random_binary failed.");

	for (int i = 0; i < params_glwe->nn * params_glwe->k; ++i)
	{
		cr_assert(sk->values[i] == 0 || sk->values[i] == 1, "Non-binary values");
	}

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Ensures glwe_sk_random_binary returns a non-NULL pointer when values != NULL.
 */
PvdaParamTest(glwe_sk_random_ternary, values_not_null, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = new_glwe_secret_key(params_glwe);

	// Draw uniformly in Zn[X] the bivGLWE secret key's values
	int status = glwe_sk_random_ternary(module, sk);

	// Asserts glwe_sk_random_uniform worked
	cr_assert(eq(int, status, 0), "glwe_sk_random_ternary failed.");

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Ensures glwe_sk_random_binary's generated key is actually binary
 */
PvdaParamTest(glwe_sk_random_ternary, values_ternary, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Create a bivGLWE secret key
	GLWESecretKey* sk = new_glwe_secret_key(params_glwe);

	int status = glwe_sk_random_ternary(module, sk);

	cr_assert(eq(int, status, 0), "glwe_sk_random_ternary failed.");

	for (int i = 0; i < params_glwe->nn * params_glwe->k; ++i)
	{
		cr_assert(sk->values[i] == 0 || sk->values[i] == 1 || sk->values[i] == -1, "Non-binary values");
	}

	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GLWE;
}

// =============================================
// |                                           |
// |     Extraction, preparation, edge cases   |
// |                                           |
// =============================================

#define KN 64
#define KK 3

static GLWEParams* key_params(void) { return new_glwe_params(KN, KK, 10, 4 * KK + 4, 0, NOISE_UNIFORM_POWER_OF_TWO); }

Test(glwe_key, new_key_is_zero_and_sized_by_the_parameters)
{
	GLWEParams* params_glwe = key_params();
	GLWESecretKey* sk       = new_glwe_secret_key(params_glwe);
	cr_assert(eq(u64, sk->nn, KN));
	cr_assert(eq(u64, sk->k, KK));
	for (int i = 0; i < KN * KK; ++i) cr_assert(eq(i64, sk->values[i], 0));
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key(NULL);
	delete_glwe_params(params_glwe);
}

Test(glwe_key, extract_poly_views_the_key)
{
	GLWEParams* params_glwe = key_params();
	GLWESecretKey* sk       = new_glwe_secret_key(params_glwe);
	for (int i = 0; i < KN * KK; ++i) sk->values[i] = i;
	for (uint64_t j = 0; j < KK; ++j)
	{
		PolyUnivZnx* s_j = glwe_sk_extract_poly(sk, j);
		cr_assert(eq(ptr, s_j, sk->values + j * KN));
		cr_assert(eq(i64, s_j[0], (int64_t)(j * KN)));
	}
	// A view: writing through it modifies the key
	glwe_sk_extract_poly(sk, KK - 1)[KN - 1] = -5;
	cr_assert(eq(i64, sk->values[KN * KK - 1], -5));
	delete_glwe_secret_key(sk);
	delete_glwe_params(params_glwe);
}

Test(glwe_key, uniform_key_bounds_and_invalid_bits)
{
	set_log_flag(0);
	GLWEParams* params_glwe = key_params();
	PvdaBackend* module     = pvda_new_spqlios_backend(KN);
	GLWESecretKey* sk       = new_glwe_secret_key(params_glwe);
	cr_assert(eq(int, glwe_sk_random_uniform(module, sk, 1), 0));  // {-1, 0}
	for (int i = 0; i < KN * KK; ++i) cr_assert(sk->values[i] == 0 || sk->values[i] == -1);
	cr_assert(eq(int, glwe_sk_random_uniform(module, sk, 64), 0));
	cr_assert(eq(int, glwe_sk_random_uniform(module, sk, 0), -1));
	cr_assert(eq(int, glwe_sk_random_uniform(module, sk, 65), -1));
	delete_glwe_secret_key(sk);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
}

Test(glwe_key, sampling_is_sized_by_the_key)
{
	// A key smaller than the ring of the module: only its own k * nn coefficients are written
	GLWEParams* params_glwe = new_glwe_params(16, 1, 10, 4, 0, NOISE_UNIFORM_POWER_OF_TWO);
	PvdaBackend* module     = pvda_new_spqlios_backend(1024);
	GLWESecretKey* sk       = new_glwe_secret_key(params_glwe);
	cr_assert(eq(int, glwe_sk_random_ternary(module, sk), 0));
	cr_assert(eq(int, glwe_sk_random_binary(module, sk), 0));
	cr_assert(eq(int, glwe_sk_random_uniform(module, sk, 3), 0));
	delete_glwe_secret_key(sk);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
}

Test(glwe_key, prepare_copies_coefficients_and_prepares_each_component)
{
	GLWEParams* params_glwe        = key_params();
	PvdaBackend* module            = pvda_new_spqlios_backend(KN);
	GLWESecretKey* sk              = new_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = new_glwe_secret_key_prepared(module, params_glwe);
	cr_assert_not_null(sk_prep);
	cr_assert(eq(ptr, (void*)sk_prep->module, module));
	cr_assert(eq(u64, sk_prep->nn, KN));
	cr_assert(eq(u64, sk_prep->k, KK));

	cr_assert(eq(int, glwe_sk_random_ternary(module, sk), 0));
	cr_assert(eq(int, glwe_sk_prepare(module, sk_prep, sk), 0));

	int64_t one[KN] = {1}, got[KN];  // the constant polynomial 1
	PolyBiv one_biv = biv_view(KN, 1, KN, one);
	PolyBiv got_biv = biv_view(KN, 1, KN, got);
	BE_PolyBivDFT* dft = pvda_new_biv_dft(module, 1);
	for (uint64_t j = 0; j < KK; ++j)
	{
		// Coefficient copy, as a view into the prepared key
		PolyUnivZnx* c = glwe_sk_prepared_extract_poly_coefs(sk_prep, j);
		cr_assert(eq(ptr, c, sk_prep->values_coef + j * KN));
		for (int p = 0; p < KN; ++p) cr_assert(eq(i64, c[p], sk->values[j * KN + p]));

		// Prepared component j: s_j * 1 = s_j
		BE_PolyUnivSvpPrep* s_j = glwe_prepared_sk_extract_poly_prep(sk_prep, j);
		cr_assert_not_null(s_j);
		cr_assert(eq(int, pvda_svp_apply_dft(module, dft, 1, s_j, &one_biv), 0));
		cr_assert(eq(int, pvda_biv_idft(module, &got_biv, dft, 1), 0));
		for (int p = 0; p < KN; ++p) cr_assert(eq(i64, got[p], sk->values[j * KN + p]), "component %lu coef %d", j, p);
	}

	pvda_delete_biv_dft(module, dft);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_glwe_secret_key_prepared(NULL);
	delete_glwe_secret_key(sk);
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
}
