#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "bivariate_polynomial.h"
#include "core/ggsw/ggsw.h"
#include "core/glwe/glwe.h"
#include "ggsw_ciphertext.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

PvdaTstParams params = {1024, 1, 4, 4, 4, -3};

/**
 * @brief Tests ggsw_secret_encrpyt
 *
 */
Test(ggsw_secret_encrypt, works)
{
	INIT_PVDA_PARAMS_GGSW(&params);

	// The decryption of a bivGLWE(m) should give m_dec = m + err, and |m_dec - m| <= 2^(-kappa*l) + 3*sigma with a 99%
	// chance
	double max_err_length      = ldexp(1.0, -params_glwe->l * params_glwe->kappa) + 3 * sigma;
	double critical_err_length = ldexp(1.0, -params_glwe->l * params_glwe->kappa) + 5 * sigma;

	GLWESecretKey* sk        = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft = alloc_glwe_secret_key_dft(params_glwe);
	GGSWCiphertext* ggsw     = new_ggsw(params_ggsw);
	PolyUniv* m_univ         = new_univ(params_glwe);
	PolyUnivDFT* m_univ_dft  = new_univ_dft(module);
	// Variables to compute the phase of each ggsw's row
	PolyBiv* phase_computed              = new_biv_poly(params_glwe);
	PolyUnivRnX* phase_observed_univ_rnx = new_univ_rnx(params_glwe);
	PolyUnivRnX* phase_expected_univ_rnx = new_univ_rnx(params_glwe);
	PolyUnivDFT* m_skj_univ_dft          = new_univ_dft(module);
	PolyUniv* m_skj_univ                 = new_univ(params_glwe);

	// Draws uniformly in (Zn[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);
	uniform_random_pol_znx(m_univ, params_glwe->nn, params_glwe->kappa);
	univ_coefs_to_dft(module, m_univ_dft, m_univ);

	ggsw_secret_encrypt(module, ggsw, sk_dft, m_univ);

	// Asserts the j-th row in the i-th PartialGGSW(m) in the ggsw is :
	// - a bivGLWE(-m * sk_j / 2^{kappa_tilde * i})), for j < k
	// - a bivGLWE(m / 2^{kappa_tilde * i}))        , for j = k
	for (uint64_t i = 1; i <= params_ggsw->l_tilde; i++)
	{
		for (uint64_t j = 0; j <= params_ggsw->k_tilde; j++)
		{
			memset(phase_computed, 0, poly_biv_bytes(params_glwe));
			memset(phase_observed_univ_rnx, 0, poly_univ_bytes(params_glwe));
			memset(m_skj_univ_dft, 0, poly_univ_bytes(params_glwe));
			memset(m_skj_univ, 0, poly_univ_bytes(params_glwe));
			memset(phase_expected_univ_rnx, 0, poly_univ_bytes(params_glwe));

			// Retrieves the phase, which should equal -m * sk_j / 2^{kappa_tilde * i}) + err
			GLWECiphertext glwe_ct = {params_glwe, ggsw_retrieve_bivglwe(ggsw, j, i)};
			glwe_secret_decrypt(module, phase_computed, sk_dft, &glwe_ct);
			biv_to_univ_rnx(params_glwe, phase_observed_univ_rnx, phase_computed);

			// Computes -m * sk_j / 2^{i*kappa_tilde}
			if (j < params_ggsw->k_tilde)
			{
				mult_vec_znx_dft(module, m_skj_univ_dft, 1, glwe_sk_extract_poly_dft(sk_dft, j), 1, m_univ_dft, 1);
				for (uint64_t p = 0; p < params_glwe->nn; p++) m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];
				pvda_vec_znx_idft(module, m_skj_univ, 1, m_skj_univ_dft, 1);
			}

			for (uint64_t p = 0; p < params_glwe->nn; p++)
				phase_expected_univ_rnx[p] =
				    ldexp((j == params_ggsw->k_tilde) ? (double)m_univ[p] : (double)m_skj_univ[p],
				          -(params_ggsw->kappa_tilde * i));

			pvda_assert_polynomial_distance(params_glwe, phase_observed_univ_rnx, phase_expected_univ_rnx,
			                                max_err_length, critical_err_length);
		}
	}

	// Clean up
	delete_univ_rnx(phase_observed_univ_rnx);
	delete_univ_dft(m_skj_univ_dft);
	delete_univ(m_skj_univ);
	delete_univ(phase_computed);
	delete_univ_rnx(phase_expected_univ_rnx);
	delete_univ(m_univ);
	delete_univ_dft(m_univ_dft);
	delete_ggsw(ggsw);
	delete_glwe_secret_key_dft(sk_dft);

	DELETE_PVDA_PARAMS_GGSW;
}
