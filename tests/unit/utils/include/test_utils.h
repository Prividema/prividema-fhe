#ifndef PVDA_UTEST_UTILS
#define PVDA_UTEST_UTILS

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
	double sigma;
} PvdaTstParams;

double generate_sigma(PvdaTstParams* p);

#define INIT_PVDA_PARAMS_GLWE(PRS)                             \
	MODULE* module          = pvda_new_module_info((PRS)->nn); \
	double sigma            = generate_sigma((PRS));           \
	GLWEParams* params_glwe = new_glwe_params((PRS)->nn, (PRS)->k, (PRS)->kappa, (PRS)->l * ((PRS)->k + 1), sigma);

#define DELETE_PVDA_PARAMS_GLWE      \
	pvda_delete_module_info(module); \
	delete_glwe_params(params_glwe);

#define INIT_PVDA_PARAMS_GGSW(PRS) \
	INIT_PVDA_PARAMS_GLWE((PRS))   \
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, (PRS)->k, (PRS)->kappa, (PRS)->l_tilde * ((PRS)->k + 1))

#define DELETE_PVDA_PARAMS_GGSW \
	DELETE_PVDA_PARAMS_GLWE     \
	delete_ggsw_params(params_ggsw);

#endif
