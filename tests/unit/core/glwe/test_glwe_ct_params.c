#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "core/glwe/glwe_ct_params.h"

#define NBASE      1024
#define KBASE      8
#define KAPPABASE  4
#define NLIMBSBASE 45
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -12

Test(new_glwe_ct_params, basic)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	cr_assert(eq(int, params != NULL, 1));
	cr_assert(eq(i64, params->N, NBASE));
	cr_assert(eq(i64, params->k, KBASE));
	cr_assert(eq(i64, params->kappa, KAPPABASE));
	cr_assert(eq(i64, params->n_limbs, NLIMBSBASE));
	cr_assert(eq(dbl, params->sigma, ldexp(1.0, SIGMABASE)));

	delete_glwe_ct_params(params);
}
