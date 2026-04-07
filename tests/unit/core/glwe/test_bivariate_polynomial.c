#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>

#include "bivariate_polynomial.h"
#include "common/rng.h"
#include "common/spqlios_alias.h"
#include "test_utils.h"
#include "univariate_polynomial.h"
#include "utils.h"

PvdaTstParams params = {1024, 1, 19, 4, 0, -7};
//! COMMON PART (begin)

Test(poly_biv_size, basic)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	cr_assert(eq(i64, glwe_params_l(params_glwe), params_glwe->l));

	DELETE_PVDA_PARAMS_GLWE;
}

Test(poly_biv_bytes, basic)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	cr_assert(eq(i64, poly_biv_bytes(params_glwe), params_glwe->nn * params_glwe->l * sizeof(int64_t)));

	DELETE_PVDA_PARAMS_GLWE;
}

Test(biv_to_univ_rnx, runs)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	PolyBiv* pol                   = new_biv_poly(params_glwe);
	PolyUnivRnX* pol_univ_computed = new_univ_rnx(params_glwe);

	//Set the first limb for the first coefficient to 1, ie, first coefficient is set to 2^-kappa
	pol[0] = 1;

	biv_to_univ_rnx(params_glwe, pol_univ_computed, pol);

	//Checks that the first coefficient is indeed 2^-kappa
	cr_assert(eq(dbl, pol_univ_computed[0], ldexp(1.0, -params_glwe->kappa)));

	//And the rest should be 0
	for (uint64_t p = 1; p < params_glwe->nn; p++)
		cr_assert(eq(dbl, pol_univ_computed[p], 0), "pol_univ_computed[%ld] = %ld ", p, pol_univ_computed[p]);

	free(pol);
	delete_univ_rnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

Test(univ_rnx_to_biv, one_test)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	PolyUnivRnX* pol_univ = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed = new_biv_poly(params_glwe);

	// Define pol_univ(X) = 2^{-params->kappa} , i.e. in Zn[X,Y] pol(X,Y) = Y
	pol_univ[0] = ldexp(1.0, -params_glwe->kappa);

	// Compute pol_univ's base-2params->kappa normalized decomposition
	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ);

	// Asserts pol_computed = Y
	cr_assert(eq(int, pol_computed[0], 1), "pol_computed[%ld, %ld] = %ld ", 0, 1, pol_computed[0]);
	for (uint64_t i = 1; i < params_glwe->l * params_glwe->nn; i++)
		cr_assert(eq(int, pol_computed[i], 0), "pol_computed[%ld, %ld] = %ld ", i / params_glwe->nn,
		          i % params_glwe->nn, pol_computed[i]);

	delete_univ_rnx(pol_univ);
	free(pol_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

Test(univ_rnx_to_biv, basic)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	PolyUnivRnX* pol_univ = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed = new_biv_poly(params_glwe);

	// Draws normaly pol_univ in Rn[X] (!= torus)
	// TODO: use a torus element instead of rnx?
	normal_random_vec(pol_univ, params_glwe->nn, 0.0, 1e-2);

	// Computes pol_univ's base-2params->kappa normalized decomposition
	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ);

	double err_length = ldexp(1.0, -(params_glwe->l - 1) * params_glwe->kappa) + 3 * DBL_EPSILON;

	// Asserts pol_computed, in Rn[X] (with Y = 2^{-params->kappa}), is equal to pol_univ
	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		double pol_computed_p = 0;
		for (uint64_t i = 1; i <= params_glwe->l; i++)
			pol_computed_p += ldexp((double)pol_computed[(i - 1) * params_glwe->nn + p], -i * params_glwe->kappa);

		//TODO: double check this
		cr_assert(epsilon_eq(dbl, torus_distance(pol_computed_p, pol_univ[p]), 0, err_length));
	}

	delete_univ_rnx(pol_univ);
	free(pol_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

Test(univ_rnx_to_biv, maths_test)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	PolyUnivRnX* pol_univ          = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed          = new_biv_poly(params_glwe);
	PolyUnivRnX* pol_univ_computed = new_univ_rnx(params_glwe);

	normal_random_vec(pol_univ, params_glwe->nn, 0.0, 1e-2);

	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ);

	biv_to_univ_rnx(params_glwe, pol_univ_computed, pol_computed);

	double err_length = ldexp(1.0, -(params_glwe->l - 1) * params_glwe->kappa) + 3 * DBL_EPSILON;

	// Asserts pol_univ_computed(X) = pol_univ(X)
	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		cr_assert(epsilon_eq(dbl, torus_distance(pol_univ[p], pol_univ_computed[p]), 0, err_length));
	}

	delete_univ_rnx(pol_univ);
	free(pol_computed);
	delete_univ_rnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

Test(univ_tnx_to_biv, maths_test)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	PolyUnivTnX* pol_univ          = new_univ_tnx(params_glwe);
	PolyBiv* pol_computed          = new_biv_poly(params_glwe);
	PolyUnivTnX* pol_univ_computed = new_univ_tnx(params_glwe);

	uniform_random_vec(params_glwe->nn, (PolyUniv*)pol_univ, 1, params_glwe->nn, 64);

	univ_tnx_to_biv(params_glwe, pol_computed, pol_univ);

	biv_to_univ_tnx(params_glwe, pol_univ_computed, pol_computed);

	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		//TODO: adjust to allow L*KAPPA < 64
		cr_assert(eq(u64, pol_univ[p], pol_univ_computed[p]));
	}

	delete_univ_tnx(pol_univ);
	free(pol_computed);
	delete_univ_tnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

//! BIV POLY PART (begin)

/**
 * @brief Tests if it returns the right size for N in [1,100]
 * @note n_limbs = (k + 1) * l
 */
Test(poly_biv_coef_number, classic_params)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	// Asserts poly_biv_coef_number returns params_glwe->nn * params_glwe->l
	cr_assert(eq(i64, poly_biv_coef_number(params_glwe), params_glwe->nn * params_glwe->l));

	DELETE_PVDA_PARAMS_GLWE;
}

Test(normal_random_biv_poly, does_not_crash)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	PolyBiv* pol = new_biv_poly(params_glwe);

	int status = normal_random_biv_poly(params_glwe, pol);

	cr_assert(eq(int, status, 0, "normal_random_biv failed."));

	free(pol);

	DELETE_PVDA_PARAMS_GLWE;
}

Test(normal_random_biv_poly, output_is_normalized)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	PolyBiv* pol = new_biv_poly(params_glwe);

	// Draw normaly pol in Zn[X,Y]
	normal_random_biv_poly(params_glwe, pol);

	// Asserts pol is normalized.
	// i.e. that each coefficient is between -2^(params->kappa-1) (inclusive) and 2^(params->kappa-1) (exculsive)
	for (uint64_t i = 0; i < params_glwe->l * params_glwe->nn; i++)
	{
		cr_assert(lt(i64, pol[i], (1 << (params_glwe->kappa - 1))));
		cr_assert(ge(i64, pol[i], -(1 << (params_glwe->kappa - 1))));
	}

	free(pol);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Test add_biv_poly correctness with random normal polynomials
 */
Test(add_biv_poly, basic)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	PolyBiv* pol_lhs      = new_biv_poly(params_glwe);
	PolyBiv* pol_rhs      = new_biv_poly(params_glwe);
	PolyBiv* sum_observed = new_biv_poly(params_glwe);

	normal_random_biv_poly(params_glwe, pol_lhs);
	normal_random_biv_poly(params_glwe, pol_rhs);

	add_biv_poly(module, params_glwe, sum_observed, pol_lhs, pol_rhs);

	// Asserts sum_computed = pol_lhs + pol_rhs
	for (uint64_t i = 0; i < params_glwe->l * params_glwe->nn; i++)
	{
		cr_assert(eq(i64, sum_observed[i], pol_lhs[i] + pol_rhs[i]));
	}

	free(pol_lhs);
	free(pol_rhs);
	free(sum_observed);

	DELETE_PVDA_PARAMS_GLWE;
}
