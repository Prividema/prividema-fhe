#include "glwegadget_utils.h"

#include <criterion/internal/new_asserts.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "glwe_ciphertext.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

void check_glwegadget(const PVDA_MODULE* module, const GLWEGadgetCiphertext* glwegad,
                      const GLWESecretKeyPrepared* sk_prep, const PolyUniv* expected, double max_err,
                      double critical_err)
{
	const GLWEGadgetParams* params_glwegad = glwegad->params;
	const GLWEParams* params_glwe          = params_glwegad->params_glwe;
	PolyBiv* phase_computed                = new_biv(params_glwe);
	PolyUnivRnX* phase_observed_univ_rnx   = new_univ_rnx(params_glwe);
	PolyUnivRnX* phase_expected_univ_rnx   = new_univ_rnx(params_glwe);

	for (uint64_t prec_lvl = 1; prec_lvl < params_glwegad->l_tilde; ++prec_lvl)
	{
		memset(phase_computed->ptr, 0, poly_biv_bytes(params_glwe));
		memset(phase_observed_univ_rnx, 0, poly_univ_bytes(params_glwe));
		memset(phase_expected_univ_rnx, 0, poly_univ_bytes(params_glwe));

		GLWECiphertext glwe_ct = {params_glwe, glwegadget_extract_bivglwe(glwegad, prec_lvl)};
		int code               = glwe_secret_decrypt(module, phase_computed, sk_prep, &glwe_ct);
		cr_assert(code == 0);
		biv_to_univ_rnx(params_glwe, phase_observed_univ_rnx, phase_computed);

		for (uint64_t p = 0; p < params_glwe->nn; p++)
			phase_expected_univ_rnx[p] = ldexp((double)expected[p], -(params_glwegad->kappa_tilde * prec_lvl));

		pvda_assert_polynomial_distance(params_glwe, phase_observed_univ_rnx, phase_expected_univ_rnx, max_err,
		                                critical_err);
	}

	delete_biv(phase_computed);
	delete_univ_rnx(phase_expected_univ_rnx);
	delete_univ_rnx(phase_observed_univ_rnx);
}

int64_t info_bits_half_prod(const GLWEParams* params_glwe, const GLWEGadgetParams* params_glwegadget)
{
	int decomp_info_bits = params_glwe->kappa * glwe_params_l_b(params_glwe);

	int log2n   = ceil(log2(params_glwe->nn));
	int log2nlt = ceil(log2(params_glwe->nn * params_glwegadget->l_tilde));

	decomp_info_bits -= log2n;
	decomp_info_bits -= log2nlt;
	decomp_info_bits -= 2;

	return decomp_info_bits;
}
