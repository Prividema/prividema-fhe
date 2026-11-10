#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <criterion/redirect.h>

#include "backend_arithmetic.h"
#include "bivariate_polynomial.h"

#include "glwe_params.h"
#include "random_polynomial.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

PvdaTstParams params = {1024, 1, 4, 4, 0, -7};
//! COMMON PART (begin)

/**
 * @brief Test biv_bytes
 */
PvdaParamTest(univ_znx_bytes, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Asserts univ_znx_bytes returns params_glwe->nn * sizeof(int64_t)
	cr_assert(eq(i64, univ_znx_bytes(params_glwe), params_glwe->nn * sizeof(int64_t)));

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Test biv_bytes
 */
PvdaParamTest(coef_dft_back_forth, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivZnx* a            = new_univ_znx(params_glwe);
	PolyUnivZnx* a_t          = new_univ_znx(params_glwe);
	BE_PolyUnivZnxDFT* res_dft = new_univ_znx_dft(module);

	univ_znx_random_uniform_pow2(module, a, params_glwe->nn, 8);

	univ_znx_to_dft(module, res_dft, a);
	univ_dft_to_znx(module, a_t, res_dft);

	for (int i = 0; i < params_glwe->nn; ++i)
	{
		cr_assert(eq(i64, a_t[i], a[i]));
	}
	delete_univ_znx(a);
	delete_univ_znx(a_t);
	delete_univ_znx_dft(module, res_dft);

	DELETE_PVDA_PARAMS_GLWE;
}

Test(tnx_rnx_encoding, known_bounded_values)
{
	PolyUnivRnX rnx_values[4] = {0, -0.25, 0.25, -0.5};
	PolyUnivTnX tnx_values[4] = {0, 0xc000000000000000L, 0x4000000000000000L, 0x8000000000000000L};

	PolyUnivRnX rnx_computed[4];
	PolyUnivTnX tnx_computed[4];

	GLWEParams* params_glwe = new_glwe_params(4, 1, 4, 2, 0.25, NOISE_UNIFORM_POWER_OF_TWO);

	univ_rnx_to_tnx(params_glwe, tnx_computed, rnx_values);
	univ_tnx_to_rnx(params_glwe, rnx_computed, tnx_values);

	for (size_t i = 0; i < sizeof(rnx_computed) / sizeof(PolyUnivRnX); ++i)
	{
		cr_assert(eq(u64, tnx_values[i], tnx_computed[i]));
		cr_assert(eq(dbl, rnx_values[i], rnx_computed[i]));
	}

	delete_glwe_params(params_glwe);
}

Test(tnx_rnx_encoding, known_outbounded_values)
{
	PolyUnivRnX rnx_values[5] = {2, -1.25, 68.25, 0.5, 0.875};
	PolyUnivTnX tnx_values[5] = {0, 0xc000000000000000L, 0x4000000000000000L, 0x8000000000000000L, 0xe000000000000000L};

	PolyUnivTnX tnx_computed[5];

	GLWEParams* params_glwe = new_glwe_params(5, 1, 4, 1, 0.25, NOISE_UNIFORM_POWER_OF_TWO);

	univ_rnx_to_tnx(params_glwe, tnx_computed, rnx_values);

	for (size_t i = 0; i < sizeof(tnx_computed) / sizeof(PolyUnivTnX); ++i)
	{
		cr_assert(eq(u64, tnx_values[i], tnx_computed[i]));
	}

	delete_glwe_params(params_glwe);
}

PvdaParamTest(tnx_rnx_encoding, back_and_forth_rnx, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);
	int vec_size = params_glwe->nn;

	PolyUnivRnX* rnx_values = new_univ_rnx(params_glwe);
	PolyUnivTnX* tnx_values = new_univ_tnx(params_glwe);
	PolyUnivRnX* rnx_final  = new_univ_rnx(params_glwe);

	rnx_random_vec(module, rnx_values, params_glwe);
	for (int i = 0; i < vec_size; ++i)
	{
		rnx_values[i] = rnx_values[i] - floor(rnx_values[i]);
	}
	univ_rnx_to_tnx(params_glwe, tnx_values, rnx_values);
	univ_tnx_to_rnx(params_glwe, rnx_final, tnx_values);

	pvda_assert_polynomial_distance(params_glwe, rnx_final, rnx_values, 3 * DBL_EPSILON, 4 * DBL_EPSILON);

	delete_univ_rnx(rnx_values);
	delete_univ_tnx(tnx_values);
	delete_univ_rnx(rnx_final);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(tnx_rnx_encoding, back_and_forth_tnx, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);
	int vec_size = params_glwe->nn;

	PolyUnivRnX* rnx_values = new_univ_rnx(params_glwe);
	PolyUnivTnX* tnx_values = new_univ_tnx(params_glwe);
	PolyUnivTnX* tnx_final  = new_univ_tnx(params_glwe);

	uint64_t precision = (1l << (64 - 53));
	univ_znx_random_uniform_pow2(module, tnx_values, vec_size, 64);

	univ_tnx_to_rnx(params_glwe, rnx_values, tnx_values);
	univ_rnx_to_tnx(params_glwe, tnx_final, rnx_values);

	for (uint64_t i = 0; i < vec_size; ++i)
	{
		uint64_t diff = tnx_final[i] > tnx_values[i] ? tnx_final[i] - tnx_values[i] : tnx_values[i] - tnx_final[i];
		cr_assert(le(u64, diff, precision));
	}

	delete_univ_rnx(rnx_values);
	delete_univ_tnx(tnx_values);
	delete_univ_tnx(tnx_final);

	DELETE_PVDA_PARAMS_GLWE;
}

// =============================================
// |                                           |
// |            Allocation and sizes           |
// |                                           |
// =============================================

Test(univ_alloc, every_type_and_size)
{
	static const uint64_t sizes[] = {1, 4, 8, 1024};
	for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s)
	{
		GLWEParams* params_glwe = new_glwe_params(sizes[s], 1, 4, 3, 0, NOISE_UNIFORM_POWER_OF_TWO);
		cr_assert(eq(u64, univ_znx_bytes(params_glwe), sizes[s] * sizeof(int64_t)));
		cr_assert(eq(u64, univ_rnx_bytes(params_glwe), sizes[s] * sizeof(double)));
		cr_assert(eq(u64, univ_tnx_bytes(params_glwe), sizes[s] * sizeof(uint64_t)));

		PolyUnivZnx* z = new_univ_znx(params_glwe);
		PolyUnivRnX* r = new_univ_rnx(params_glwe);
		PolyUnivTnX* t = new_univ_tnx(params_glwe);
		cr_assert_not_null(z, "nn = %lu", sizes[s]);
		cr_assert_not_null(r);
		cr_assert_not_null(t);
		// ZnX polynomials are 64-byte aligned for the backends, even when nn * 8 < 64
		cr_assert(eq(u64, (uintptr_t)z % 64, 0));
		// The whole polynomial is writable
		for (uint64_t i = 0; i < sizes[s]; ++i) z[i] = r[i] = t[i] = 0;

		delete_univ_znx(z);
		delete_univ_rnx(r);
		delete_univ_tnx(t);
		delete_glwe_params(params_glwe);
	}
	// Deleting NULL does nothing
	delete_univ_znx(NULL);
	delete_univ_rnx(NULL);
	delete_univ_tnx(NULL);
}

// =============================================
// |                                           |
// |                 DFT / SVP                 |
// |                                           |
// =============================================

/** @brief Naive negacyclic product in Z[X]/(X^n + 1). */
static void naive_negacyclic(int64_t* res, const int64_t* a, const int64_t* b, uint64_t n)
{
	memset(res, 0, n * sizeof(int64_t));
	for (uint64_t i = 0; i < n; ++i)
		for (uint64_t j = 0; j < n; ++j)
		{
			int64_t prod = a[i] * b[j];
			if (i + j < n) res[i + j] += prod;
			else res[i + j - n] -= prod;
		}
}

Test(univ_dft, round_trip_edge_values)
{
	const uint64_t nn       = 64;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 4, 3, 0, NOISE_UNIFORM_POWER_OF_TWO);
	PolyUnivZnx* a          = new_univ_znx(params_glwe);
	PolyUnivZnx* back       = new_univ_znx(params_glwe);
	BE_PolyUnivZnxDFT* dft  = new_univ_znx_dft(module);
	cr_assert_not_null(dft);

	// Zero, a single monomial, and large alternating values (the DFT is exact for |coefficients| well below 2^50)
	for (int c = 0; c < 3; ++c)
	{
		for (uint64_t i = 0; i < nn; ++i)
			a[i] = c == 0 ? 0 : (c == 1 ? (i == nn - 1 ? -1 : 0) : ((i & 1) ? -(INT64_C(1) << 40) : (INT64_C(1) << 40) - 1));
		cr_assert(eq(int, univ_znx_to_dft(module, dft, a), 0));
		cr_assert(eq(int, univ_dft_to_znx(module, back, dft), 0));
		for (uint64_t i = 0; i < nn; ++i) cr_assert(eq(i64, back[i], a[i]), "case %d coef %lu", c, i);
	}

	delete_univ_znx_dft(module, dft);
	delete_univ_znx(back);
	delete_univ_znx(a);
	delete_glwe_params(params_glwe);
	pvda_delete_backend(module);
}

Test(univ_svp_prep, prepared_product_is_negacyclic_product)
{
	const uint64_t nn       = 64;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 4, 3, 0, NOISE_UNIFORM_POWER_OF_TWO);
	PolyUnivZnx* s          = new_univ_znx(params_glwe);
	int64_t a[64], got[64], expected[64];
	PolyBiv a_biv           = biv_view(nn, 1, nn, a);
	PolyBiv got_biv         = biv_view(nn, 1, nn, got);
	BE_PolyUnivSvpPrep* s_p = new_univ_svp_prep(module);
	BE_PolyBivDFT* res_dft  = pvda_new_biv_dft(module, 1);
	cr_assert_not_null(s_p);

	cr_assert(univ_znx_random_ternary(module, s, nn) == 0);
	cr_assert(univ_znx_random_uniform_pow2(module, a, nn, 12) == 0);
	cr_assert(eq(int, univ_znx_to_svp_prep(module, s_p, s), 0));
	cr_assert(eq(int, pvda_svp_apply_dft(module, res_dft, 1, s_p, &a_biv), 0));
	cr_assert(eq(int, pvda_biv_idft(module, &got_biv, res_dft, 1), 0));

	naive_negacyclic(expected, s, a, nn);
	for (uint64_t i = 0; i < nn; ++i) cr_assert(eq(i64, got[i], expected[i]), "coef %lu", i);

	pvda_delete_biv_dft(module, res_dft);
	delete_univ_svp_prep(module, s_p);
	delete_univ_znx(s);
	delete_glwe_params(params_glwe);
	pvda_delete_backend(module);
}

// =============================================
// |                                           |
// |          TnX <-> RnX edge cases           |
// |                                           |
// =============================================

Test(tnx_rnx_encoding, rnx_edge_values)
{
	GLWEParams* params_glwe = new_glwe_params(6, 1, 4, 3, 0, NOISE_UNIFORM_POWER_OF_TWO);
	// Tiny negative values reduce to 1 - epsilon, which rounds to 1.0 in double: they must give 0, not overflow
	PolyUnivRnX rnx[6]       = {-1e-30, -0x1p-80, 1e-30, 0.5, -0.5, 0x1p-64};
	const PolyUnivTnX exp[6] = {0, 0, 0, UINT64_C(1) << 63, UINT64_C(1) << 63, 1};
	PolyUnivTnX tnx[6];
	cr_assert(eq(int, univ_rnx_to_tnx(params_glwe, tnx, rnx), 0));
	for (int i = 0; i < 6; ++i) cr_assert(eq(u64, tnx[i], exp[i]), "coef %d", i);
	delete_glwe_params(params_glwe);
}

Test(tnx_rnx_encoding, tnx_edge_values)
{
	GLWEParams* params_glwe = new_glwe_params(5, 1, 4, 3, 0, NOISE_UNIFORM_POWER_OF_TWO);
	// Read as signed numerators over 2^64: the result is in [-1/2, 1/2)
	PolyUnivTnX tnx[5]       = {0, 1, UINT64_MAX, UINT64_C(1) << 63, (UINT64_C(1) << 63) - 1};
	const PolyUnivRnX exp[5] = {0.0, 0x1p-64, -0x1p-64, -0.5, 0.5};  // 2^63 - 1 rounds to 1/2 in double
	PolyUnivRnX rnx[5];
	cr_assert(eq(int, univ_tnx_to_rnx(params_glwe, rnx, tnx), 0));
	for (int i = 0; i < 5; ++i) cr_assert(eq(dbl, rnx[i], exp[i]), "coef %d", i);
	delete_glwe_params(params_glwe);
}

// =============================================
// |                                           |
// |              Print functions              |
// |                                           |
// =============================================

Test(univ_print, znx, .init = cr_redirect_stdout)
{
	PolyUnivZnx pol[4] = {3, -2, 0, INT64_MIN};
	print_poly_univ_ZnX(pol, 4);
	PolyUnivZnx single[1] = {-7};
	print_poly_univ_ZnX(single, 1);
	fflush(stdout);
	cr_assert_stdout_eq_str("3 - 2 X^1 + 0 X^2 - 9223372036854775808 X^3\n-7\n");
}

Test(univ_print, rnx, .init = cr_redirect_stdout)
{
	PolyUnivRnX pol[3] = {0.5, -1.25, 2};
	print_poly_univ_RnX(pol, 3);
	fflush(stdout);
	cr_assert_stdout_eq_str("0.5 - 1.25 X^1 + 2 X^2\n");
}

Test(univ_print, tnx, .init = cr_redirect_stdout)
{
	PolyUnivTnX pol[2] = {1, UINT64_MAX};
	print_poly_univ_TnX(pol, 2);
	fflush(stdout);
	cr_assert_stdout_eq_str("1 + 18446744073709551615 X^1\n");
}

Test(univ_print, znx_mp, .init = cr_redirect_stdout)
{
	PolyUnivZnxMP* mp = new_univ_znx_mp(2, 2);
	PolyUnivZnx v[2]  = {-1, 5};
	univ_znx_to_znx_mp(mp, v);
	print_poly_univ_ZnX_mp(mp);
	fflush(stdout);
	cr_assert_stdout_eq_str("X^0: 0xffffffffffffffffffffffffffffffff\nX^1: 0x00000000000000000000000000000005\n");
	delete_univ_znx_mp(mp);
}
