#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "core/glwe/glwe_params.h"

#define NBASE      1024
#define KBASE      8
#define KAPPABASE  4
#define NLIMBSBASE 45
#define LBASE      NLIMBSBASE / (KBASE + 1)

Test(new_glwe_ct_params, power_of_two_noise)
{
	GLWEParams* params = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, 3 - LBASE * KAPPABASE),
	                                     NOISE_UNIFORM_POWER_OF_TWO);

	cr_assert(eq(int, params != NULL, 1));
	cr_assert(eq(i64, params->nn, NBASE));
	cr_assert(eq(i64, params->k, KBASE));
	cr_assert(eq(i64, params->kappa, KAPPABASE));
	cr_assert(eq(i64, params->ciphertext_nb_limbs, NLIMBSBASE));

	//For NOISE_UNIFORM_POWER_OF_TWO, we should have ceil(log2(sqrt(3), sigma*2^-(K*l)))
	//In other words, the log2 rounded up of sqrt(3) times the standard deviation rescaled to the last limb
	cr_assert(eq(dbl, params->fast_uniform_nb_bits, 4));

	delete_glwe_params(params);
}
