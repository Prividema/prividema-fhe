#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>

#include "common/rng.h"
#include "common/spqlios_alias.h"
#include "glwe_params.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

PvdaTstParams params = {1024, 1, 4, 4, 0, -7};
//! COMMON PART (begin)

/**
 * @brief Test poly_biv_bytes
 */
PvdaParamTest(poly_univ_bytes, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Asserts poly_univ_bytes returns params_glwe->nn * sizeof(int64_t)
	cr_assert(eq(i64, poly_univ_bytes(params_glwe), params_glwe->nn * sizeof(int64_t)));

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Test poly_biv_bytes
 */
PvdaParamTest(coef_dft_back_forth, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUniv* a          = new_univ(params_glwe);
	PolyUniv* a_t        = new_univ(params_glwe);
	PolyUnivDFT* res_dft = new_univ_dft(module);

	uniform_random_vec(params_glwe->nn, a, 1, params_glwe->nn, 8);

	univ_coefs_to_dft(module, res_dft, a);
	univ_dft_to_coefs(module, a_t, res_dft);

	for (int i = 0; i < params_glwe->nn; ++i)
	{
		cr_assert(eq(i64, a_t[i], a[i]));
	}
	delete_univ(a);
	delete_univ(a_t);
	delete_univ_dft(res_dft);

	DELETE_PVDA_PARAMS_GLWE;
}

Test(tnx_rnx_encoding, known_bounded_values)
{
	PolyUnivRnX rnx_values[4] = {0, -0.25, 0.25, -0.5};
	PolyUnivTnX tnx_values[4] = {0, 0xc000000000000000L, 0x4000000000000000L, 0x8000000000000000L};

	PolyUnivRnX rnx_computed[4];
	PolyUnivTnX tnx_computed[4];

	GLWEParams* params_glwe = new_glwe_params(4, 1, 4, 1, 0);

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

	GLWEParams* params_glwe = new_glwe_params(5, 1, 4, 1, 0);

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

	normal_random_vec(rnx_values, vec_size, 0, 0.2);
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
	uniform_random_pol_znx(tnx_values, vec_size, 64);

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
