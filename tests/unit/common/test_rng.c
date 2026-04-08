#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdlib.h>

#include "common/logger.h"
#include "common/rng.h"
#include "core/glwe/glwe_params.h"
#include "stat_utils.h"
#include "test_utils.h"

#define NB_SAMPLES 1000

// -------------------------------------------------------------------------------------
// Uniform Distribution Test
//
// We use the Chi-Squared test : https://en.wikipedia.org/wiki/Chi-squared_test
// -------------------------------------------------------------------------------------

/**
 * Function that apply a Chi-squared test.
 *
 * It splits the interval [-2^nb_bits, 2^nb_bits) into
 * a given number of boxes and then build a distribution
 * by generating random numbers.
 */
int rand_uniform_aux(uint64_t nb_bits, int nb_boxes)
{
	int64_t max, min;
	if (nb_bits == 8 * sizeof(int64_t))
	{
		max = INT64_MAX;
		min = INT64_MIN;
	}
	else
	{
		max = (1 << (nb_bits - 1)) - 1;
		min = -(1 << (nb_bits - 1));
	}

	int64_t boxes[nb_boxes];
	memset(boxes, 0, nb_boxes * sizeof(boxes[0]));

	double step = ((double)max - (double)min) / (double)nb_boxes;
	for (size_t i = 0; i < NB_SAMPLES; i++)
	{
		int64_t sample = 0;
		if (rand_uniform(&sample, nb_bits) < 0) return -1;

		// Epsilon 1e9 used to avoid floating point precision bad rounding
		int box_n = (int)(((double)sample - (double)min + 1e-9) / step);
		if (box_n == nb_boxes) --box_n;

		++boxes[box_n];
	}

	// Apply Chi squared test
	double expected = (double)NB_SAMPLES / (double)nb_boxes;  // For a uniform distribution
	double tt       = 0.0;
	for (int i = 0; i < nb_boxes; i++)
	{
		double num = (((double)boxes[i]) - expected);
		tt += (num * num) / expected;
	}

	// The test has 5% chance of failing
	if (tt < chi_critical_05[nb_boxes - 2]) return 0;
	return 1;
}

// Test rand_uniform on the whole interval

Test(rand_uniform_64, test_rand_uniform)
{
	int count = 0;
	for (int i = 0; i < 100; i++)
	{
		int c;
		if ((c = rand_uniform_aux(64, 100)) < 0) cr_fatal("Error occured during generation");
		count += c;
	}

	// As Chi-squared has 5% chance of failing, we count the number of times it fails.
	// On 100 iterations it fails 5 times in average.
	// We determined experimentally the standard deviation sigma = 4
	// So the number of failures should be in range 5 +- 4.
	cr_assert(ge(int, count, 1), "The number of errors should be between in range 5 +- 4");
	cr_assert(le(int, count, 9), "The number of errors should be between in range 5 +- 4");
}

// Test rand_uniform on the interval [-32768,32767]
Test(rand_uniform_16, test_rand_uniform)
{
	int count = 0;
	for (int i = 0; i < 100; i++) count += rand_uniform_aux(16, 100);

	// As Chi-squared has 5% chance of failing, we count the number of times it fails.
	// On 100 iterations it fails 5 times in average.
	// We determined experimentally the standard deviation sigma = 4
	// So the number of failures should be in range 5 +- 4.
	cr_assert(ge(int, count, 1), "The number of errors should be between in range 5 +- 4");
	cr_assert(le(int, count, 9), "The number of errors should be between in range 5 +- 4");
}

// Test rand_normal with Jarque-Bera test
Test(rand_normal, test_rand_normal)
{
	double data[NB_SAMPLES];
	for (int i = 0; i < NB_SAMPLES; i++)
		if (rand_normal(data + i, 0, 1) < 0) cr_fail("rand_normal failed");

	double jjbb = jarque_bera(data, NB_SAMPLES);
	cr_assert(lt(dbl, jjbb, chi_critical_05[1]), "Expect %f < %f\n", jjbb, chi_critical_05[1]);
}

PvdaTstParams params = {4, 2, 4, 1, 1, -1};

PvdaParamTest(normal_random_vec, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	//TODO: RnX vec of size params_glwe->nn*params_glwe->k??
	VecUnivRnX* pol_univ = malloc(params_glwe->nn * params_glwe->k * sizeof(double));
	cr_assert(pol_univ != NULL);

	cr_assert(normal_random_vec(pol_univ, params_glwe->nn * params_glwe->k, 0.0, 0.001) == 0);

	free(pol_univ);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(uniform_random_vec_dft, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);
	VecUnivDFT* res_dft = malloc(params_glwe->nn * params_glwe->k * sizeof(double));

	VecUniv* res = malloc(params_glwe->nn * params_glwe->k * sizeof(int64_t));

	cr_assert(res != NULL);
	cr_assert(res_dft != NULL);

	cr_assert(uniform_random_vec_znx_dft(module, res_dft, params_glwe->k, 2) == 0);

	pvda_vec_znx_dft(module, res_dft, params_glwe->k, res, params_glwe->k, params_glwe->nn);

	free(res_dft);
	free(res);

	DELETE_PVDA_PARAMS_GLWE;
}
