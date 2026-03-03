#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "core/glwe/glwe_ciphertext.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "utils.h"
#include "vec_znx_arithmetic_private.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  4
#define NLIMBSBASE (KBASE + 1) * 2
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  1e-7

//! COMMON PART (begin)

/**
 * @brief Tests whether glwe_size computes the right size of a GLWE ciphertext.
 */
Test(glwe_size, basic)
{
	GLWECtParams* params_glwe= new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

	cr_assert(eq(i64, glwe_size(params_glwe), NLIMBSBASE));

	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Tests whether glwe_bytes computes the right number of bytes in a GLWE ciphertext.
 */
Test(glwe_bytes, basic)
{
	GLWECtParams* params_glwe= new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

	cr_assert(eq(i64, glwe_bytes(params_glwe), NLIMBSBASE * NBASE * 8));

	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Tests whether mult_vec_znx_dft multiply correctly two Zn[X] polynomials a and b. Ie res = a*b.
 */
Test(mult_vec_znx_dft, size_equal_one)
{
	GLWECtParams* params_glwe= new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module       = new_module_info(NBASE, FFT64);

	int64_t* res         = calloc(poly_univ_bytes(params_glwe), 1);

	// def a = 1 + X
	int64_t* a = calloc(poly_univ_bytes(params_glwe), 1);
	inplace_uniform_random_vec(NBASE, a, 1, NBASE, 14);

	// def b = 1 + X
	int64_t* b = calloc(poly_univ_bytes(params_glwe), 1);
	inplace_uniform_random_vec(NBASE, b, 1, NBASE, 14);

	double* res_dft = calloc(poly_univ_bytes(params_glwe), 1);
	double* a_dft   = calloc(poly_univ_bytes(params_glwe), 1);
	double* b_dft   = calloc(poly_univ_bytes(params_glwe), 1);

	vec_znx_dft_p(module, res_dft, 1, res, 1, NBASE);
	vec_znx_dft_p(module, a_dft, 1, a, 1, NBASE);
	vec_znx_dft_p(module, b_dft, 1, b, 1, NBASE);

	// res_dft = DFT(a*b)
	mult_vec_znx_dft(module, res_dft, 1, a_dft, 1, b_dft, 1);

	// res = a*b
	vec_znx_idft_p(module, res, 1, res_dft, 1);

	// Compare the real coefficient res_p for each p in [0, NBASE -1] with the res_p mult_vec_znx_dft computed
	// coefficient.
	for (int64_t p = 0; p < NBASE; p++) {
		int64_t acc = 0;
		for (uint64_t k = 0; k <= p; k++) {
			acc += a[k] * b[p - k];
		}
		for (uint64_t k = p + 1; k < NBASE; k++) {
			acc += -a[k] * b[NBASE + p - k];
		}
		cr_log_info("acc %ld res[p] %ld X^%ld", acc, res[p], p);
		cr_assert(eq(i64, res[p], acc));
	}

	free(res);
	free(res_dft);
	free(a);
	free(a_dft);
	free(b);
	free(b_dft);
	delete_module_info(module);
	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Tests whether mult_vec_znx_dft correctly multiplies two Zn[X] vectors a and b component-wise..
 * It draws a random uniform size, ie a random uniform number of Zn[X] polynomials.
 * Ie a = (a_i), b = (b_i) -> res = (a_i * b_i). Where a_i and b_i are in Zn[X]
 */
Test(mult_vec_znx_dft, random_size)
{
	GLWECtParams* params_glwe= new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module       = new_module_info(NBASE, FFT64);

	int64_t size         = 0;

	while (size <= 0) {
		rand_uniform(&size, 8);
	}

	int64_t* res = calloc(poly_univ_bytes(params_glwe) * size, 1);

	int64_t* a   = calloc(poly_univ_bytes(params_glwe) * size, 1);
	inplace_uniform_random_vec(NBASE, a, size, NBASE, 14);

	int64_t* b = calloc(poly_univ_bytes(params_glwe) * size, 1);
	inplace_uniform_random_vec(NBASE, b, size, NBASE, 14);

	double* res_dft = calloc(poly_univ_bytes(params_glwe) * size, 1);
	double* a_dft   = calloc(poly_univ_bytes(params_glwe) * size, 1);
	double* b_dft   = calloc(poly_univ_bytes(params_glwe) * size, 1);

	vec_znx_dft_p(module, res_dft, size, res, size, NBASE);
	vec_znx_dft_p(module, a_dft, size, a, size, NBASE);
	vec_znx_dft_p(module, b_dft, size, b, size, NBASE);

	// res_dft = DFT(a*b)
	mult_vec_znx_dft(module, res_dft, size, a_dft, size, b_dft, size);

	// res = a*b
	vec_znx_idft_p(module, res, size, res_dft, size);

	// Compare the real coefficient res_p for each p in [0, NBASE -1] with the res_p mult_vec_znx_dft computed
	// coefficient.
	for (int64_t i = 0; i < size; i++) {
		for (int64_t p = 0; p < NBASE; p++) {
			int64_t acc = 0;
			for (uint64_t k = 0; k <= p; k++) {
				acc += a[i * NBASE + k] * b[i * NBASE + p - k];
			}
			for (uint64_t k = p + 1; k < NBASE; k++) {
				acc += -a[i * NBASE + k] * b[i * NBASE + NBASE + p - k];
			}
			cr_log_info("acc %ld res[i*NBASE + p] %ld X^%ld", acc, res[i * NBASE + p], p);
			cr_assert(eq(i64, res[i * NBASE + p], acc));
		}
	}

	free(res);
	free(res_dft);
	free(a);
	free(a_dft);
	free(b);
	free(b_dft);
	delete_module_info(module);
	delete_glwe_ct_params(params_glwe);
}

//! GLWE PART (begin)

/**
 * @brief Tests whether glwe_bytes computes the right number of coefficient in a GLWE ciphertext.
 */
Test(glwe_coef_number, basic)
{
	GLWECtParams* params_glwe= new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

	cr_assert(eq(i64, glwe_coef_number(params_glwe), NLIMBSBASE * NBASE));

	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Tests whether new_glwe returns a non-NULL pointer.
 */
Test(new_glwe, basic)
{
	GLWECtParams* params_glwe= new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GLWECiphertext* ct   = new_glwe(params_glwe);

	cr_assert(eq(int, (ct != NULL) && (ct->vec != NULL), 1));

	delete_glwe(ct);
	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Tests whether add_glwe adds two GLWE ciphertexts.
 */
Test(add_glwe, basic)
{
	GLWECtParams* params_glwe= new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

	GLWECiphertext* ct_l = new_glwe(params_glwe);
	GLWECiphertext* ct_r = new_glwe(params_glwe);
	GLWECiphertext* res  = new_glwe(params_glwe);

	inplace_uniform_random_vec(NBASE, ct_l->vec, params_glwe->n_limbs, NBASE, KAPPABASE - 1);
	inplace_uniform_random_vec(NBASE, ct_r->vec, params_glwe->n_limbs, NBASE, KAPPABASE - 1);

	add_glwe(res, ct_l, ct_r);

	for (int64_t i = 1; i < LBASE; i++)
		for (int64_t j = 0; j < KBASE + 1; j++)
			for (int64_t p = 0; p < NBASE; p++)
				cr_assert(eq(res->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p],
				             ct_l->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p] +
				                 ct_r->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p]));

	delete_glwe(ct_l);
	delete_glwe(ct_r);
	delete_glwe(res);
	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Tests whether const_mult_glwe multiply a GLWE ciphertext by a ZnX polynomial.
 */
Test(const_mult_glwe, without_normalization)
{
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module       = new_module_info(NBASE, FFT64);

	GLWECiphertext* res  = new_glwe(params_glwe);

	// Draws uniformly the GLWE ciphertext and the ZnX polynomial
	GLWECiphertext* ct = new_glwe(params_glwe);
	inplace_uniform_random_vec(NBASE, ct->vec, params_glwe->n_limbs, NBASE, KAPPABASE - 1);

	PolyUniv* u        = new_uniform_random_vec(NBASE, KAPPABASE - 1);
	PolyUnivDFT* u_dft = malloc(NBASE * sizeof(int64_t));
	vec_znx_dft_p(module, u_dft, 1, u, 1, NBASE);

	const_mult_glwe(module, res, u_dft, ct, 0);

	for (int64_t i = 1; i <= LBASE; i++)
		for (int64_t j = 0; j < KBASE + 1; j++) {
			PolyUniv* ct_ij = ct->vec + (i - 1) * (KBASE + 1) * NBASE + j * NBASE;
			for (int64_t p = 0; p < NBASE; p++) {
				int64_t acc = 0;
				for (uint64_t k = 0; k <= p; k++) {
					acc += u[k] * ct_ij[p - k];
				}
				for (uint64_t k = p + 1; k < NBASE; k++) {
					acc += -u[k] * ct_ij[NBASE + p - k];
				}
				cr_assert(eq(i64, res->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p], acc));
			}
		}

	free(u);
	free(u_dft);
	delete_module_info(module);
	delete_glwe(ct);
	delete_glwe(res);
	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Tests whether const_mult_glwe multiply a GLWE ciphertext by a ZnX polynomial.
 */
Test(const_mult_glwe, with_normalization)
{
	GLWECtParams* params_glwe= new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module       = new_module_info(NBASE, FFT64);

	GLWECiphertext* res  = new_glwe(params_glwe);

	// Draws uniformly the GLWE ciphertext and the ZnX polynomial
	GLWECiphertext* ct = new_glwe(params_glwe);
	inplace_uniform_random_vec(NBASE, ct->vec, params_glwe->n_limbs, NBASE, KAPPABASE - 1);

	PolyUniv* u        = new_uniform_random_vec(NBASE, KAPPABASE - 1);
	PolyUnivDFT* u_dft = malloc(NBASE * sizeof(int64_t));
	vec_znx_dft_p(module, u_dft, 1, u, 1, NBASE);

	const_mult_glwe(module, res, u_dft, ct, 1);

	for (int64_t j = 0; j < KBASE + 1; j++)
		for (int64_t p = 0; p < NBASE; p++)
			for (int64_t i = 1; i <= LBASE; i++) {
				int64_t remainder = 0;
				for (int64_t i = LBASE; i >= 1; i--) {
					PolyUniv* ct_ij = ct->vec + (i - 1) * (KBASE + 1) * NBASE + j * NBASE;

					int64_t acc     = 0;
					for (uint64_t k = 0; k <= p; k++) {
						acc += u[k] * ct_ij[p - k];
					}
					for (uint64_t k = p + 1; k < NBASE; k++) {
						acc += -u[k] * ct_ij[NBASE + p - k];
					}

					cr_assert(
					    eq(i64,
					       (res->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p] - (acc + remainder)) % KAPPABASE,
					       0),
					    "Equality failed at j = %ld p = %ld i = %ld with acc = %ld reminder = %ld and res = %ld", j, p,
					    i, acc, remainder, res->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p]);

					remainder = acc >= 0 ? (acc + (1 << KAPPABASE - 1)) / (1 << KAPPABASE)
					                     : (acc - (1 << KAPPABASE - 1) + 1) / (1 << KAPPABASE);
				}
			}

	free(u);
	free(u_dft);
	delete_module_info(module);
	delete_glwe(ct);
	delete_glwe(res);
	delete_glwe_ct_params(params_glwe);
}

//! GLWE IN DFT PART (begin)

/**
 * @brief Tests whether glwe_coef_number_dft computes the right number of coefficient in a GLWE ciphertext in DFT space.
 */
Test(glwe_coef_number_dft, basic)
{
	GLWECtParams* params_glwe= new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

	cr_assert(eq(i64, glwe_coef_number_dft(params_glwe), NLIMBSBASE * NBASE / 2));

	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Tests whether new_glwe_dft returns a non-NULL pointer.
 */
Test(new_glwe_dft, basic)
{
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GLWECiphertextDFT* ct = new_glwe_dft(params_glwe);

	cr_assert(eq(int, (ct != NULL) && (ct->vec != NULL), 1));

	delete_glwe_dft(ct);
	delete_glwe_ct_params(params_glwe);
}

Test(add_glwe_dft, basic)
{
	GLWECtParams* params_glwe         = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module                = new_module_info(NBASE, FFT64);

	GLWECiphertextDFT* ct_l_dft   = new_glwe_dft(params_glwe);
	GLWECiphertextDFT* ct_r_dft   = new_glwe_dft(params_glwe);
	GLWECiphertextDFT* ct_sum_dft = new_glwe_dft(params_glwe);

	inplace_uniform_random_vec_znx_dft(module, ct_l_dft->vec, params_glwe->n_limbs, KAPPABASE - 1);
	inplace_uniform_random_vec_znx_dft(module, ct_r_dft->vec, params_glwe->n_limbs, KAPPABASE - 1);

	add_glwe_dft(ct_sum_dft, ct_l_dft, ct_r_dft);

	for (int64_t i = 1; i < LBASE; i++)
		for (int64_t j = 0; j < KBASE + 1; j++)
			for (int64_t p = 0; p < NBASE; p++)
				cr_assert(eq(ct_sum_dft->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p],
				             ct_l_dft->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p] +
				                 ct_r_dft->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p]));

	delete_module_info(module);
	delete_glwe_dft(ct_l_dft);
	delete_glwe_dft(ct_r_dft);
	delete_glwe_dft(ct_sum_dft);
	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Tests whether const_mult_glwe_dft multiply a GLWE ciphertext by a ZnX polynomial.
 */
Test(const_mult_glwe_dft, without_normalization)
{
	GLWECtParams* params_glwe      = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module             = new_module_info(NBASE, FFT64);

	GLWECiphertextDFT* res_dft = new_glwe_dft(params_glwe);
	VecBiv* res_vec            = malloc(glwe_bytes(params_glwe));

	// Draws uniformly the GLWE ciphertext and computes it out of DFT space
	GLWECiphertextDFT* ct_dft = new_glwe_dft(params_glwe);
	inplace_uniform_random_vec_znx_dft(module, ct_dft->vec, params_glwe->n_limbs, KAPPABASE - 1);

	VecBiv* ct_vec = malloc(glwe_bytes(params_glwe));
	vec_znx_idft_p(module, ct_vec, glwe_size(params_glwe), ct_dft->vec, glwe_size(params_glwe));

	// Draws uniformly the ZnX polynomial and computes it ouf of DFT space
	PolyUniv* u        = new_uniform_random_vec(NBASE, KAPPABASE - 1);
	PolyUnivDFT* u_dft = malloc(NBASE * sizeof(int64_t));
	vec_znx_dft_p(module, u_dft, 1, u, 1, NBASE);

	const_mult_glwe_dft(module, res_dft, u_dft, ct_dft, 0);

	// Computes res out of DFT space
	vec_znx_idft_p(module, res_vec, glwe_size(params_glwe), res_dft->vec, glwe_size(params_glwe));

	for (int64_t i = 1; i <= LBASE; i++)
		for (int64_t j = 0; j < KBASE + 1; j++) {
			PolyUniv* ct_ij = ct_vec + (i - 1) * (KBASE + 1) * NBASE + j * NBASE;
			for (int64_t p = 0; p < NBASE; p++) {
				int64_t acc = 0;
				for (uint64_t k = 0; k <= p; k++) {
					acc += u[k] * ct_ij[p - k];
				}
				for (uint64_t k = p + 1; k < NBASE; k++) {
					acc += -u[k] * ct_ij[NBASE + p - k];
				}
				cr_assert(eq(i64, res_vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p], acc));
			}
		}

	free(u);
	free(u_dft);
	free(ct_vec);
	free(res_vec);
	delete_module_info(module);
	delete_glwe_dft(ct_dft);
	delete_glwe_dft(res_dft);
	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Tests whether const_mult_glwe_dft multiply a GLWE ciphertext by a ZnX polynomial.
 */
Test(const_mult_glwe_dft, with_normalization)
{
	GLWECtParams* params_glwe      = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module             = new_module_info(NBASE, FFT64);

	GLWECiphertextDFT* res_dft = new_glwe_dft(params_glwe);
	VecBiv* res_vec            = malloc(glwe_bytes(params_glwe));

	// Draws uniformly the GLWE ciphertext and computes it out of DFT space
	GLWECiphertextDFT* ct_dft = new_glwe_dft(params_glwe);
	inplace_uniform_random_vec_znx_dft(module, ct_dft->vec, params_glwe->n_limbs, KAPPABASE - 1);

	VecBiv* ct_vec = malloc(glwe_bytes(params_glwe));
	vec_znx_idft_p(module, ct_vec, glwe_size(params_glwe), ct_dft->vec, glwe_size(params_glwe));

	// Draws uniformly the ZnX polynomial and computes it ouf of DFT space
	PolyUniv* u        = new_uniform_random_vec(NBASE, KAPPABASE - 1);
	PolyUnivDFT* u_dft = malloc(NBASE * sizeof(int64_t));
	vec_znx_dft_p(module, u_dft, 1, u, 1, NBASE);

	const_mult_glwe_dft(module, res_dft, u_dft, ct_dft, 1);

	// Computes res out of DFT space
	vec_znx_idft_p(module, res_vec, glwe_size(params_glwe), res_dft->vec, glwe_size(params_glwe));

	for (int64_t j = 0; j < KBASE + 1; j++)
		for (int64_t p = 0; p < NBASE; p++)
			for (int64_t i = 1; i <= LBASE; i++) {
				int64_t remainder = 0;
				for (int64_t i = LBASE; i >= 1; i--) {
					PolyUniv* ct_ij = ct_vec + (i - 1) * (KBASE + 1) * NBASE + j * NBASE;

					int64_t acc     = 0;
					for (uint64_t k = 0; k <= p; k++) {
						acc += u[k] * ct_ij[p - k];
					}
					for (uint64_t k = p + 1; k < NBASE; k++) {
						acc += -u[k] * ct_ij[NBASE + p - k];
					}

					cr_assert(
					    eq(i64,
					       (res_vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p] - (acc + remainder)) % KAPPABASE, 0),
					    "Equality failed at j = %ld p = %ld i = %ld with acc = %ld reminder = %ld and res = %ld", j, p,
					    i, acc, remainder, res_vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p]);

					remainder = acc >= 0 ? (acc + (1 << KAPPABASE - 1)) / (1 << KAPPABASE)
					                     : (acc - (1 << KAPPABASE - 1) + 1) / (1 << KAPPABASE);
				}
			}

	free(u);
	free(u_dft);
	free(ct_vec);
	free(res_vec);
	delete_module_info(module);
	delete_glwe_dft(ct_dft);
	delete_glwe_dft(res_dft);
	delete_glwe_ct_params(params_glwe);
}