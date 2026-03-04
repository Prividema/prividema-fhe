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
	PolyUniv** values = new_glwe_secret_key_values(NBASE, KBASE);

	cr_assert(eq(int, values != NULL, 1));
	for (uint64_t j = 0; j < KBASE; j++) cr_assert(eq(int, values[j] != NULL, 1));

	delete_glwe_secret_key_values(values, KBASE);
}

/**
 * @brief Ensure new_uniform_secret_key_values creates no NULL-pointer.
 */
Test(new_uniform_glwe_secret_key_values, basic)
{
	MODULE* module    = new_module_info(NBASE, FFT64);
	PolyUniv** values = new_uniform_glwe_secret_key_values(module, KBASE, 3);

	cr_assert(eq(int, values != NULL, 1));
	for (uint64_t j = 0; j < KBASE; j++) cr_assert(eq(int, values[j] != NULL, 1));

	delete_module_info(module);
	delete_glwe_secret_key_values(values, KBASE);
}

/**
 * @brief Tests wether transform_secret_key_values_dft_to_not_dft actually transforms the values of the secret key in
 * DFT space, out of DFT space.
 */
Test(transform_glwe_secret_key_values_dft_to_not_dft, basic)
{
	MODULE* module           = new_module_info(NBASE, FFT64);
	PolyUnivDFT** values_dft = new_uniform_glwe_secret_key_values_dft(module, KBASE, 3);
	PolyUniv** values        = transform_glwe_secret_key_values_dft_to_not_dft(module, (const PolyUnivDFT**)values_dft, KBASE);

	for (uint64_t j = 0; j < KBASE; j++)
		for (uint64_t p = 0; p < NBASE; p++) cr_log_info("%" PRId64 " X^%" PRId64, values[j][p], p);

	delete_module_info(module);
	delete_glwe_secret_key_values_dft(values_dft, KBASE);
	delete_glwe_secret_key_values(values, KBASE);
}

/**
 * @brief Ensures new_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
Test(new_glwe_secret_key, values_not_null)
{
	GLWESecretKey* sk = new_glwe_secret_key(NBASE, KBASE);

	cr_assert(eq(int, sk != NULL, 1));

	cr_assert(eq(int, sk->values != NULL, 1));

	for (uint64_t j = 0; j < KBASE; j++) cr_assert(eq(int, sk->values[j] != NULL, 1));

	delete_glwe_secret_key(sk);
}

/**
 * @brief Ensures new_uniform_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
Test(new_uniform_glwe_secret_key, values_not_null)
{
	MODULE* module    = new_module_info(NBASE, FFT64);
	GLWESecretKey* sk = new_uniform_glwe_secret_key(module, KBASE, 2);

	cr_assert(eq(int, sk != NULL, 1));
	cr_assert(eq(int, sk->values != NULL, 1));

	for (uint64_t j = 0; j < KBASE; j++) cr_assert(eq(int, sk->values[j] != NULL, 1));

	delete_module_info(module);
	delete_glwe_secret_key(sk);
}

/**
 * @brief Tests wether transform_glwe_secret_key_dft_to_not_dft transforms the secret key in the DFT domain, out of DFT
 * space.
 */
Test(transform_glwe_secret_key_dft_to_not_dft, basic)
{
	MODULE* module           = new_module_info(NBASE, FFT64);
	GLWESecretKeyDFT* sk_dft = new_uniform_glwe_secret_key_dft(module, KBASE, 3);
	GLWESecretKey* sk        = transform_glwe_secret_key_dft_to_not_dft(module, sk_dft);

	for (uint64_t j = 0; j < KBASE; j++)
		for (uint64_t p = 0; p < NBASE; p++) cr_log_info("%" PRId64 " X^%" PRId64, sk->values[j][p], p);

	delete_module_info(module);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_secret_key(sk);
}

//! GLWE KEY PART IN DFT SPACE (begin)

/**
 * @brief Ensure new_secret_key_values_dft creates no NULL-pointer.
 */
Test(new_glwe_secret_key_values_dft, basic)
{
	PolyUnivDFT** values_dft = new_glwe_secret_key_values_dft(NBASE, KBASE);

	cr_assert(eq(int, values_dft != NULL, 1));
	for (uint64_t j = 0; j < KBASE; j++) {
		cr_assert(eq(int, values_dft[j] != NULL, 1));
	}

	delete_glwe_secret_key_values_dft(values_dft, KBASE);
}

/**
 * @brief Ensure new_uniform_secret_key_values_dft creates no NULL-pointer.
 */
Test(new_uniform_glwe_secret_key_values_dft, basic)
{
	MODULE* module           = new_module_info(NBASE, FFT64);
	PolyUnivDFT** values_dft = new_uniform_glwe_secret_key_values_dft(module, KBASE, 3);

	cr_assert(eq(int, values_dft != NULL, 1));
	for (uint64_t j = 0; j < KBASE; j++) cr_assert(eq(int, values_dft[j] != NULL, 1));

	delete_module_info(module);
	delete_glwe_secret_key_values_dft(values_dft, KBASE);
}

/**
 * @brief Tests wether transform_secret_key_values_not_dft_to_dft actually transforms the values of the secret key out
 * of DFT space, in the DFT domain.
 */
Test(transform_glwe_secret_key_values_not_dft_to_dft, basic)
{
	MODULE* module     = new_module_info(NBASE, FFT64);
	PolyUniv** values        = new_uniform_glwe_secret_key_values(module, KBASE, 3);
	PolyUnivDFT** values_dft = transform_glwe_secret_key_values_not_dft_to_dft(module, (const PolyUniv**)values, KBASE);

	for (uint64_t j = 0; j < KBASE; j++)
		for (uint64_t p = 0; p < NBASE; p++) cr_log_info("%lf X^%" PRId64, values_dft[j][p], p);

	delete_module_info(module);
	delete_glwe_secret_key_values(values, KBASE);
	delete_glwe_secret_key_values_dft(values_dft, KBASE);
}

/**
 * @brief Ensures new_glwe_secret_key_dft returns a non-NULL pointer when values != NULL.
 */
Test(new_glwe_secret_key_dft, values_not_null)
{
	GLWESecretKeyDFT* sk_dft = new_glwe_secret_key_dft(NBASE, KBASE);

	cr_assert(eq(int, sk_dft != NULL, 1));
	cr_assert(eq(int, sk_dft->values != NULL, 1));

	for (uint64_t j = 0; j < KBASE; j++) {
		cr_assert(eq(int, sk_dft->values[j] != NULL, 1));
	}

	delete_glwe_secret_key_dft(sk_dft);
}

/**
 * @brief Tests whether new_uniform_glwe_secret_key_gen works as intended.
 */
Test(new_uniform_glwe_secret_key_dft, what_s_inside)
{
	MODULE* module           = new_module_info(NBASE, FFT64);

	GLWESecretKeyDFT* sk_dft = new_uniform_glwe_secret_key_dft(module, KBASE, 2);

	cr_assert(eq(int, sk_dft != NULL, 1));
	cr_assert(eq(int, sk_dft->values != NULL, 1));

	for (uint64_t j = 0; j < KBASE; j++) cr_assert(eq(int, sk_dft->values[j] != NULL, 1));

	delete_module_info(module);
	delete_glwe_secret_key_dft(sk_dft);
}

/**
 * @brief Tests wether transform_glwe_secret_key_not_dft_to_dft transforms the secret key out of DFT space, in DFT
 * space.
 */
Test(transform_glwe_secret_key_not_dft_to_dft, basic)
{
	MODULE* module           = new_module_info(NBASE, FFT64);
	GLWESecretKey* sk        = new_uniform_glwe_secret_key(module, KBASE, 3);
	GLWESecretKeyDFT* sk_dft = transform_glwe_secret_key_not_dft_to_dft(module, sk);

	for (uint64_t j = 0; j < KBASE; j++)
		for (uint64_t p = 0; p < NBASE; p++) cr_log_info("%lf X^%" PRId64, sk_dft->values[j][p], p);

	delete_module_info(module);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
}
