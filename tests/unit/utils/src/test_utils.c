
#include "test_utils.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <criterion/parameterized.h>
#include <math.h>

#include "glwe_params.h"
#include "math.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

double generate_sigma(PvdaTstParams* p)
{
	if (p->sigma > 0.0) return p->sigma;
	if (p->sigma < 0.0) return ldexp(1.0, (int)p->sigma);

	// If p == 0 defaut sigma which has arbitrarily chosen to be 4 bits of the last limb: 2^(4-K*l_a)
	int64_t l_a = ((p->ciphertext_nb_limbs + 1) / (p->k + 1));
	return ldexp(1.0, 4 - l_a * p->kappa);
}

void pvda_assert_polynomial_distance(const GLWEParams* params_glwe, PolyUnivRnX* a, PolyUnivRnX* b, double max_err,
                                     double critical_err)
{
	const int prob_factor = 3;

	int big_error_count = 0;
	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		double diff = rnx_torus_distance(a[p], b[p]);

		cr_assert(lt(dbl, diff, critical_err), "Difference outside range (too big)");

		int cond = diff < max_err;

		if (!cond) big_error_count++;
	}

	int max_fails = (int)(prob_factor * 0.0027 * (double)params_glwe->nn);
	cr_assert(big_error_count <= max_fails, "Too many values not following the dist");
}

void assert_tnx_close_enough(uint64_t a, uint64_t b, uint64_t bits)
{
	uint64_t diff = tnx_torus_distance(a, b);
	if (bits >= 64)
	{
		uint64_t max_diff = 1;
		cr_assert(le(u64, diff, max_diff));
	}
	else
	{
		uint64_t max_diff = 1ULL << (64 - bits);
		cr_assert(lt(u64, diff, max_diff));
	}
}

struct criterion_test_params default_params_fn()
{
	static PvdaTstParams default_params[] = {
	    {.nn                        = 1024,
	     .k                         = 1,
	     .kappa                     = 4,
	     .ciphertext_nb_limbs       = 8l * 2,
	     .ciphertext_nb_limbs_tilde = 8l * 2,
	     .sigma                     = 0},  // toy params, let default sigma
	    {.nn                        = (1 << 14),
	     .k                         = 1,
	     .kappa                     = 19,
	     .ciphertext_nb_limbs       = 15l * 2,
	     .ciphertext_nb_limbs_tilde = 15l * 2,
	     .sigma                     = 0},  // lattigo params, default sigma
	    {.nn                        = 1024,
	     .k                         = 4,
	     .kappa                     = 8,
	     .ciphertext_nb_limbs       = 9l * 5,
	     .ciphertext_nb_limbs_tilde = 9l * 5,
	     .sigma                     = 0},  // k > 1 params
	    {.nn                        = 1024,
	     .k                         = 4,
	     .kappa                     = 8,
	     .ciphertext_nb_limbs       = 9l * 5 - 1,
	     .ciphertext_nb_limbs_tilde = 9l * 5 - 1,
	     .sigma                     = 0},  // k > 1 l_a != l_b params
	    {.nn                        = 1024,
	     .k                         = 4,
	     .kappa                     = 8,
	     .ciphertext_nb_limbs       = 9l * 5 - 1,
	     .ciphertext_nb_limbs_tilde = 9l * 5,
	     .sigma                     = 0},  // k > 1 l_a != l_b params

	};

	return cr_make_param_array(PvdaTstParams, default_params, sizeof(default_params) / sizeof(default_params[0]));
}

int rnx_random_vec(PolyUnivRnX* res, GLWEParams* params_glwe)
{
	int status = -1;

	PolyUnivTnX* tmp_tnx = new_univ_tnx(params_glwe);

	// (ab)use the fact that tnx and Z mod 2^64 are isomorphic and
	// the memory representation is the same for isomprphic values
	uniform_random_pol_znx((PolyUniv*)tmp_tnx, params_glwe->nn, 64);

	univ_tnx_to_rnx(params_glwe, res, tmp_tnx);

	status = 0;
cleanup:
	free(tmp_tnx);
	return status;
}
/*

OLD CODE that might be useful later on

void printf_glwe(GLWECiphertext* glwe)
{
    for (int64_t j = 0; j < KBASE + 1; j++) printf_poly_biv(glwe->vec + j * NBASE, (KBASE + 1) * NBASE, NBASE, LBASE);
}

void printf_glwe_dft(MODULE* module, GLWECiphertextDFT* glwe_dft)
{
    VecBiv* glwe_vec = malloc(glwe_params_bytes(glwe_dft->params));
    pvda_vec_znx_idft(module, glwe_vec, glwe_params_n_limbs(glwe_dft->params), glwe_dft->vec,
                      glwe_params_n_limbs(glwe_dft->params));

    for (int64_t j = 0; j < KBASE + 1; j++) printf_poly_biv(glwe_vec + j * NBASE, (KBASE + 1) * NBASE, NBASE, LBASE);

    free(glwe_vec);
}

*/
