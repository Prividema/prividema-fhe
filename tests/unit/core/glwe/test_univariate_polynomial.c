#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <inttypes.h>
#include <math.h>

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
Test(poly_univ_bytes, basic)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	// Asserts poly_univ_bytes returns params_glwe->nn * sizeof(int64_t)
	cr_assert(eq(i64, poly_univ_bytes(params_glwe), params_glwe->nn * sizeof(int64_t)));

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Test poly_biv_bytes
 */
Test(coef_dft_back_forth, basic)
{
	INIT_PVDA_PARAMS_GLWE(&params);

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
