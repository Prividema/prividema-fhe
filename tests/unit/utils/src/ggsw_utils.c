
#include <criterion/internal/new_asserts.h>
#include <string.h>

#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe_arithmetic.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "glwegadget_ciphertext.h"
#include "spqlios_alias.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

void check_ggsw(const PvdaBackend* module, const GGSWCiphertext* ggsw, const GLWESecretKeyPrepared* sk_prep,
                const PolyUniv* expected, double max_err, double critical_err)
{
	const GGSWParams* params_ggsw        = ggsw->params;
	const GLWEParams* params_glwe        = params_ggsw->params_glwe;
	PolyBiv* phase_computed              = new_biv(params_glwe);
	PolyUnivRnX* phase_observed_univ_rnx = new_univ_rnx(params_glwe);
	PolyUnivRnX* phase_expected_univ_rnx = new_univ_rnx(params_glwe);
	PolyUnivDFT* m_skj_univ_dft          = new_univ_dft(module);
	PolyUnivDFT* m_univ_dft              = new_univ_dft(module);
	PolyUniv* m_skj_univ                 = new_univ(params_glwe);

	univ_coefs_to_dft(module, m_univ_dft, expected);

	for (uint64_t ij = 0; ij < ggsw_num_rows(params_ggsw); ++ij)
	{
		uint64_t sk_idx   = (ij % (params_ggsw->k_tilde + 1));
		uint64_t prec_lvl = ij / (params_ggsw->k_tilde + 1) + 1;
		memset(phase_computed->ptr, 0, poly_biv_bytes(params_glwe));
		memset(phase_observed_univ_rnx, 0, poly_univ_bytes(params_glwe));
		memset(m_skj_univ_dft, 0, poly_univ_bytes(params_glwe));
		memset(m_skj_univ, 0, poly_univ_bytes(params_glwe));
		memset(phase_expected_univ_rnx, 0, poly_univ_bytes(params_glwe));

		// Retrieves the phase, which should equal -m * sk_j / 2^{kappa_tilde * i}) + err
		GLWECiphertext glwe_ct = {params_glwe, ggsw_retrieve_bivglwe(ggsw, sk_idx, prec_lvl)};
		int code               = glwe_secret_decrypt(module, phase_computed, sk_prep, &glwe_ct);
		cr_assert(code == 0);
		biv_to_univ_rnx(params_glwe, phase_observed_univ_rnx, phase_computed, 0);

		// Computes -m * sk_j / 2^{i*kappa_tilde}
		if (sk_idx < params_ggsw->k_tilde)
		{
			pvda_svp_apply_dft_to_dft(module, m_skj_univ_dft, 1, glwe_prepared_sk_extract_poly_dft(sk_prep, sk_idx),
			                          m_univ_dft, 1);
			for (uint64_t p = 0; p < params_glwe->nn; p++) m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];
			PolyBiv m_skj_biv = new_biv_view(params_glwe->nn, 1, params_glwe->nn, m_skj_univ);
			pvda_vec_znx_idft(module, &m_skj_biv, m_skj_univ_dft, 1);
		}

		for (uint64_t p = 0; p < params_glwe->nn; p++)
			phase_expected_univ_rnx[p] =
			    ldexp((sk_idx == params_ggsw->k_tilde) ? (double)expected[p] : (double)m_skj_univ[p],
			          -(params_ggsw->kappa_tilde * prec_lvl));

		pvda_assert_polynomial_distance(params_glwe, phase_observed_univ_rnx, phase_expected_univ_rnx, max_err,
		                                critical_err);
	}

	delete_biv(phase_computed);
	delete_univ_rnx(phase_observed_univ_rnx);
	delete_univ_rnx(phase_expected_univ_rnx);
	delete_univ_dft(module, m_skj_univ_dft);
	delete_univ_dft(module, m_univ_dft);
	delete_univ(m_skj_univ);
}
