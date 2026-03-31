#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <inttypes.h>
#include <math.h>

#include "common/rng.h"
#include "common/spqlios_alias.h"
#include "glwe_params.h"
#include "univariate_polynomial.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  4
#define NLIMBSBASE (KBASE + 1) * 4
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -7

//! COMMON PART (begin)

/**
 * @brief Test poly_biv_bytes
 */
Test(poly_univ_bytes, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Asserts poly_univ_bytes returns NBASE * sizeof(int64_t)
	cr_assert(eq(i64, poly_univ_bytes(params_glwe), NBASE * sizeof(int64_t),
	             "poly_univ_bytes failed: got %" PRId64 ", expected %" PRId64, poly_univ_bytes(params_glwe),
	             NBASE * sizeof(int64_t)));

	// Clean up
	delete_glwe_params(params_glwe);
}

/**
 * @brief Test poly_biv_bytes
 */
Test(coef_dft_back_forth, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module          = pvda_new_module_info(NBASE);
	PolyUniv* a             = new_univ(params_glwe);
	PolyUniv* a_t           = new_univ(params_glwe);
	PolyUnivDFT* res_dft    = new_univ_dft(module);

	uniform_random_vec(NBASE, a, 1, NBASE, 8);

	univ_coefs_to_dft(module, res_dft, a);
	univ_dft_to_coefs(module, a_t, res_dft);

	for (int i = 0; i < NBASE; ++i)
	{
		cr_assert(eq(i64, a_t[i], a[i]));
	}
	// Clean up
	delete_glwe_params(params_glwe);
	delete_univ(a);
	delete_univ(a_t);
	delete_univ_dft(res_dft);
	pvda_delete_module_info(module);
}
