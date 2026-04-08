#ifndef PVDA_UTEST_UTILS
#define PVDA_UTEST_UTILS

#include <criterion/parameterized.h>
#include <stdint.h>
#include <sys/types.h>

#include "glwe_params.h"
#include "rng.h"

/*
 * Asserts (and fails a test if not met) the following condition:
 *
 * The difference between polynomials a and b coefficient-wise is over max_err at most 3*0.27% of the time.
 * The difference between polynomials a and b is, coefficient-wise, never equal or greater to max_err.
 *
 * @param params_glwe The GLWE params
 * @param a           A RnX polynomial
 * @param b           A RnX polynomial
 * @param stdev_diff  The standard deviation (sigma) of the aforementioned Normal dist.
 *
 *
 */

void pvda_assert_polynomial_distance(const GLWEParams* params_glwe, PolyUnivRnX* a, PolyUnivRnX* b, double max_err,
                                     double critical_err);

typedef struct pvda_tst_params_t
{
	uint64_t nn;
	uint64_t k;
	uint64_t kappa;
	uint64_t l;
	uint64_t l_tilde;

	/** Can be :
	 * - =0 if one wants the default computation according to the other params
	 * - <0 to specify a power of 2
	 * - >0 to specify the exact value
	 */
	double sigma;
} PvdaTstParams;

/**
 * Internal use function to fill the sigma value of PvdaTstParams and get its value.
 *
 * Since sigma can be specified as a power, the actual stdev or left to a default to be
 * computed according to the other params, this function does that.
 *
 *
 */
double generate_sigma(PvdaTstParams* p);

#define INIT_PVDA_PARAMS_BASE(PRS) MODULE* module = pvda_new_module_info((PRS)->nn);

#define DELETE_PVDA_PARAMS_BASE    pvda_delete_module_info(module);

#define INIT_PVDA_PARAMS_GLWE(PRS)                   \
	INIT_PVDA_PARAMS_BASE((PRS))                     \
	double sigma            = generate_sigma((PRS)); \
	GLWEParams* params_glwe = new_glwe_params((PRS)->nn, (PRS)->k, (PRS)->kappa, (PRS)->l * ((PRS)->k + 1), sigma);

#define DELETE_PVDA_PARAMS_GLWE \
	DELETE_PVDA_PARAMS_BASE     \
	delete_glwe_params(params_glwe);

#define INIT_PVDA_PARAMS_GGSW(PRS) \
	INIT_PVDA_PARAMS_GLWE((PRS))   \
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, (PRS)->k, (PRS)->kappa, (PRS)->l_tilde * ((PRS)->k + 1))

#define DELETE_PVDA_PARAMS_GGSW \
	DELETE_PVDA_PARAMS_GLWE     \
	delete_ggsw_params(params_ggsw);

/**
 * Sample default parameters generator function for Prividema parametrized tests.
 *
 * This function uses hard-coded parameter tuples (structs), but other functions
 * might want to generate the parameters at execution time (if doing cartesian products, for example)
 *
 */
struct criterion_test_params default_params_fn();

/**
 * Macro to define a parametrized test for Prividema.
 *
 * It always uses the PvdaTstParams object to pass the parameters
 *
 * @param suite_name   The suite name as in normal Criterion
 * @param test_name    The test name as in normal Criterion
 * @param generator_fn The name of a function returning a struct criterion_test_params with parameters to use for the test
 *                     See default_params_fn for an example of a static generator function.
 *
 *
 *
 */
#define PvdaParamTest(suite_name, test_name, generator_fn)                        \
	ParameterizedTestParameters(suite_name, test_name) { return generator_fn(); } \
                                                                                  \
	ParameterizedTest(PvdaTstParams* param, suite_name, test_name)

#endif
