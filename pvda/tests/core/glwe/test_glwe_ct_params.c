#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <math.h>

#include "logger.h"

#include "core/glwe/glwe_params.h"

#define NBASE      1024
#define KBASE      8
#define KAPPABASE  4
#define NLIMBSBASE 45
#define LBASE      NLIMBSBASE / (KBASE + 1)

Test(new_glwe_ct_params, power_of_two_noise)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, 3 - LBASE * KAPPABASE),
	                                          NOISE_UNIFORM_POWER_OF_TWO);

	cr_assert(eq(int, params_glwe != NULL, 1));
	cr_assert(eq(i64, params_glwe->nn, NBASE));
	cr_assert(eq(i64, params_glwe->k, KBASE));
	cr_assert(eq(i64, params_glwe->kappa, KAPPABASE));
	cr_assert(eq(i64, params_glwe->ciphertext_nb_limbs, NLIMBSBASE));

	//For NOISE_UNIFORM_POWER_OF_TWO, we should have ceil(log2(sqrt(3), sigma*2^-(K*l)))
	//In other words, the log2 rounded up of sqrt(3) times the standard deviation rescaled to the last limb
	cr_assert(eq(dbl, params_glwe->fast_uniform_nb_bits, 4));

	delete_glwe_params(params_glwe);
}

// =============================================
// |                                           |
// |          Layout and derived sizes         |
// |                                           |
// =============================================

Test(glwe_params_layout, l_a_l_b_and_sizes)
{
	// n_limbs = k * l_a + l_b, with l_b = l_a (n_limbs = 15) and l_b = l_a - 1 (n_limbs = 14), k = 2
	GLWEParams* eq_params = new_glwe_params(64, 2, 10, 15, 0, NOISE_UNIFORM_POWER_OF_TWO);
	cr_assert_not_null(eq_params);
	cr_assert(eq(u64, glwe_params_l_a(eq_params), 5));
	cr_assert(eq(u64, glwe_params_l_b(eq_params), 5));
	cr_assert(eq(u64, glwe_params_n_limbs(eq_params), 15));
	cr_assert(eq(u64, glwe_coef_number(eq_params), 15 * 64));
	cr_assert(eq(u64, glwe_coef_number_dft(eq_params), 15 * 64 / 2));
	cr_assert(eq(u64, glwe_params_bytes(eq_params), 15 * 64 * sizeof(int64_t)));
	cr_assert(eq(dbl, glwe_bivariate_epsilon(eq_params), ldexp(1.0, -50)));

	GLWEParams* short_b = new_glwe_params(64, 2, 10, 14, 0, NOISE_UNIFORM_POWER_OF_TWO);
	cr_assert_not_null(short_b);
	cr_assert(eq(u64, glwe_params_l_a(short_b), 5));
	cr_assert(eq(u64, glwe_params_l_b(short_b), 4));
	cr_assert(eq(dbl, glwe_bivariate_epsilon(short_b), ldexp(1.0, -40)));

	// Smallest valid layout: one limb for every polynomial
	GLWEParams* tiny = new_glwe_params(1, 1, 1, 2, 0, NOISE_UNIFORM_POWER_OF_TWO);
	cr_assert_not_null(tiny);
	cr_assert(eq(u64, glwe_params_l_a(tiny), 1));
	cr_assert(eq(u64, glwe_params_l_b(tiny), 1));

	delete_glwe_params(eq_params);
	delete_glwe_params(short_b);
	delete_glwe_params(tiny);
	delete_glwe_params(NULL);
}

Test(glwe_params_layout, rejects_invalid_parameters)
{
	set_log_flag(0);
	cr_assert_null(new_glwe_params(0, 1, 10, 4, 0, NOISE_UNIFORM_POWER_OF_TWO));     // nn = 0
	cr_assert_null(new_glwe_params(64, 0, 10, 4, 0, NOISE_UNIFORM_POWER_OF_TWO));    // k = 0
	cr_assert_null(new_glwe_params(64, 1, 0, 4, 0, NOISE_UNIFORM_POWER_OF_TWO));     // kappa = 0
	cr_assert_null(new_glwe_params(64, 1, 64, 4, 0, NOISE_UNIFORM_POWER_OF_TWO));    // kappa > 63
	cr_assert_null(new_glwe_params(64, 1, 10, 0, 0, NOISE_UNIFORM_POWER_OF_TWO));    // no limbs
	cr_assert_null(new_glwe_params(64, 2, 10, 13, 0, NOISE_UNIFORM_POWER_OF_TWO));   // l_a = 4, l_b = 5
	cr_assert_null(new_glwe_params(64, 4, 10, 42, 0, NOISE_UNIFORM_POWER_OF_TWO));   // l_a = 8, l_b = 10
	cr_assert_null(new_glwe_params(64, 1, 10, 4, -1.0, NOISE_UNIFORM_POWER_OF_TWO)); // negative sigma
	cr_assert_null(new_glwe_params(64, 1, 10, 4, NAN, NOISE_NORMAL));                // NaN sigma
	cr_assert_null(new_glwe_params(64, 1, 10, 4, INFINITY, NOISE_NORMAL));
	cr_assert_null(new_glwe_params(64, 1, 10, 4, 0.0, (NoiseType)42));               // unknown noise type
	// kappa = 63 is the largest supported base
	GLWEParams* p = new_glwe_params(64, 1, 63, 2, 0, NOISE_UNIFORM_POWER_OF_TWO);
	cr_assert_not_null(p);
	delete_glwe_params(p);
}

// =============================================
// |                                           |
// |            Noise of each type             |
// |                                           |
// =============================================

// K * l_a = 40: sigma = 2^-30 is 2^10 on the last limb
#define NOISE_PARAMS(sigma, type) new_glwe_params(64, 1, 10, 8, (sigma), (type))

Test(glwe_params_noise, sigma_must_fit_on_the_last_limb)
{
	set_log_flag(0);
	cr_assert_null(NOISE_PARAMS(ldexp(1.0, -41), NOISE_NORMAL));  // < 1 on the last limb
	cr_assert_null(NOISE_PARAMS(ldexp(1.0, 24), NOISE_NORMAL));   // > 2^63 on the last limb
	GLWEParams* p = NOISE_PARAMS(ldexp(1.0, -40), NOISE_NORMAL);  // exactly 1
	cr_assert_not_null(p);
	delete_glwe_params(p);
}

Test(glwe_params_noise, uniform_power_of_two)
{
	GLWEParams* p = NOISE_PARAMS(ldexp(1.0, -30), NOISE_UNIFORM_POWER_OF_TWO);
	// ceil(log2(sqrt(3) * 2^10)) = 11
	cr_assert(eq(u64, p->fast_uniform_nb_bits, 11));
	cr_assert(isnan(glwe_params_stdev(p)));
	delete_glwe_params(p);

	p = NOISE_PARAMS(0, NOISE_UNIFORM_POWER_OF_TWO);  // no noise
	cr_assert(eq(u64, p->fast_uniform_nb_bits, 0));
	delete_glwe_params(p);
}

Test(glwe_params_noise, uniform)
{
	GLWEParams* p = NOISE_PARAMS(ldexp(1.0, -30), NOISE_UNIFORM);
	cr_assert(eq(i64, p->uniform_range, llround(sqrt(3.0) * 1024)));
	double b = (double)p->uniform_range;
	cr_assert(eq(dbl, glwe_params_stdev(p), ldexp(sqrt(b * (b + 1) / 3), -40)));
	cr_assert(fabs(glwe_params_stdev(p) / ldexp(1.0, -30) - 1) < 1e-3);
	delete_glwe_params(p);
}

Test(glwe_params_noise, binomial)
{
	set_log_flag(0);
	GLWEParams* p = NOISE_PARAMS(ldexp(1.0, -30), NOISE_BINOMIAL);
	// n = 4 sigma^2 = 2^22 trials, stdev sqrt(n) / 2 = 2^10 exactly
	cr_assert(eq(u64, p->binomial_n, UINT64_C(1) << 22));
	cr_assert(eq(dbl, glwe_params_stdev(p), ldexp(1.0, -30)));
	delete_glwe_params(p);

	// n is rounded up to an even number
	p = NOISE_PARAMS(ldexp(1.5, -40), NOISE_BINOMIAL);  // 4 sigma^2 = 9 -> 10
	cr_assert(eq(u64, p->binomial_n, 10));
	delete_glwe_params(p);

	// 2^30 trials is the limit
	p = NOISE_PARAMS(ldexp(1.0, -26), NOISE_BINOMIAL);  // n = 2^30
	cr_assert_not_null(p);
	delete_glwe_params(p);
	cr_assert_null(NOISE_PARAMS(ldexp(1.0, -25), NOISE_BINOMIAL));  // n = 2^32
}

Test(glwe_params_noise, normal)
{
	GLWEParams* p = NOISE_PARAMS(ldexp(1.0, -30), NOISE_NORMAL);
	cr_assert(eq(dbl, p->normal_sigma, ldexp(1.0, -30)));
	// The stdev is on the torus scale, not rescaled to the last limb
	cr_assert(eq(dbl, glwe_params_stdev(p), ldexp(1.0, -30)));
	delete_glwe_params(p);
}
