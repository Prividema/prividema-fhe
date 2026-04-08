#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdlib.h>

#include "common/logger.h"
#include "common/rng.h"
#include "core/glwe/glwe_params.h"
#include "test_utils.h"

#define NB_SAMPLES 1000

// Chi-squared critical value for alpha = 0.05 depending on the degrees of freedom
double chi_critical_05[100] = {
    3.841,   5.991,   7.815,   9.488,   11.070,  12.592,  14.067,  15.507,  16.919,  18.307,  19.675,  21.026,  22.362,
    23.685,  24.996,  26.296,  27.587,  28.869,  30.144,  31.410,  32.671,  33.924,  35.172,  36.415,  37.652,  38.885,
    40.113,  41.337,  42.557,  43.773,  44.985,  46.194,  47.400,  48.602,  49.802,  50.998,  52.192,  53.384,  54.572,
    55.758,  56.942,  58.124,  59.304,  60.481,  61.656,  62.830,  64.001,  65.171,  66.339,  67.505,  68.669,  69.832,
    70.993,  72.153,  73.311,  74.468,  75.624,  76.778,  77.931,  79.082,  80.232,  81.381,  82.529,  83.675,  84.821,
    85.965,  87.108,  88.250,  89.391,  90.531,  91.670,  92.808,  93.945,  95.081,  96.217,  97.351,  98.484,  99.617,
    100.749, 101.879, 103.010, 104.139, 105.267, 06.395,  107.522, 108.648, 109.773, 110.898, 112.022, 113.145, 114.268,
    115.390, 116.511, 117.632, 118.752, 119.871, 120.990, 122.108, 123.225, 124.342};

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

// -------------------------------------------------------------------------------------
// Gaussian Distribution Test
//
// We use the Jarque-Bera test : https://en.wikipedia.org/wiki/Jarque%E2%80%93Bera_test
// -------------------------------------------------------------------------------------

// Jarque-Bera test
double jarque_bera(const double* x, int n)
{
	if (n < 3) return NAN;

	double mean = 0.0;
	for (int i = 0; i < n; i++) mean += x[i];
	mean /= n;

	// Compute 2nd, 3rd, 4th moments
	double m2 = 0.0, m3 = 0.0, m4 = 0.0;
	for (int i = 0; i < n; i++)
	{
		double d  = x[i] - mean;
		double d2 = d * d;
		m2 += d2;
		m3 += d2 * d;
		m4 += d2 * d2;
	}
	m2 /= n;
	m3 /= n;
	m4 /= n;

	// Skewness and kurtosis
	double ss = m3 / pow(m2, 1.5);
	double kk = m4 / (m2 * m2);

	// Jarque-Bera statistic
	return (n / 6.0) * (ss * ss + ((kk - 3.0) * (kk - 3.0)) / 4.0);
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

//---------------------------------------------
//
//-----------------------------------------

PvdaTstParams params = {4, 2, 4, 1, 1, -1};

Test(normal_random_vec, basic)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	//TODO: RnX vec of size params_glwe->nn*params_glwe->k??
	VecUnivRnX* pol_univ = malloc(params_glwe->nn * params_glwe->k * sizeof(double));
	cr_assert(pol_univ != NULL);

	cr_assert(normal_random_vec(pol_univ, params_glwe->nn * params_glwe->k, 0.0, 0.001) == 0);

	free(pol_univ);

	DELETE_PVDA_PARAMS_GLWE;
}

Test(uniform_random_vec_dft, basic)
{
	INIT_PVDA_PARAMS_GLWE(&params);
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
