#include "glwe.h"

#include <sys/types.h>

#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "logger.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

//! bivGLWE PART (begin)
int glwe_secret_encrypt_phase(const MODULE* module, GLWECiphertext* glwe, const GLWESecretKeyDFT* sk_dft,
                              const PolyBiv* phase)
{
	int status = -1;

	const GLWEParams* params = (const GLWEParams*)glwe->params;
	// GLWE parameters
	uint64_t nn      = params->nn;
	uint64_t k       = params->k;
	uint64_t kappa   = params->kappa;
	uint64_t n_limbs = params->n_limbs;
	uint64_t l       = n_limbs / (k + 1);

	// The acc variable is used for the sum
	// As a reminder, if we deonte acc_j its value at cycle j,
	// at cycle j + 1 we compute
	// acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
	// In other words, acc is <A, SK>
	PolyBiv* acc         = new_biv_poly(params);
	PolyBivDFT* as_j_dft = new_biv_poly_dft(params);  // DFT(sk_j) * DFT(a_j)
	PolyBiv* as_j        = NULL;

	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_masking");
	CHECK_ALLOC(as_j_dft, "as_j_dft's calloc failed in glwe_secret_masking");

	// Draws uniformly in Zn[X,Y] the ajs'
	CHECK_CALL(uniform_random_vec(k * nn, glwe->vec, l, (k + 1) * nn, kappa),
	           "A generation failed in glwe_secret_masking_dft");

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	if (k > 1)
	{
		as_j = new_biv_poly(params);  // sk_j * a_j
		CHECK_ALLOC(as_j, "as_j's calloc failed in glwe_secret_masking");

		for (uint64_t j = 0; j < k; j++)
		{
			// The j-th component of the DFT encoding of the secret key
			PolyUnivDFT* sk_j_univ_dft = glwe_sk_extract_poly_dft(sk_dft, j);

			// Computes DFT(sk_j) * DFT(a_j)
			pvda_svp_apply_dft(module, as_j_dft, l, sk_j_univ_dft, glwe_extract_start_poly(glwe, j), l, (k + 1) * nn);

			// Undo DFT to retreive sk_j * a_j
			CHECK_CALL(pvda_vec_znx_idft(module, as_j, l, as_j_dft, l),
			           "vec_znx_idft_p failed in glwe_secret_masking_ggsw_lib");

			add_biv_poly(module, params, acc, acc, as_j);
		}
	}
	else
	{
		PolyUnivDFT* sk_j_univ_dft = glwe_sk_extract_poly_dft(sk_dft, 0);

		pvda_svp_apply_dft(module, as_j_dft, l, sk_j_univ_dft, glwe_extract_start_poly(glwe, 0), l, (k + 1) * nn);

		CHECK_CALL(pvda_vec_znx_idft(module, acc, l, as_j_dft, l),
		           "vec_znx_idft_p failed in glwe_secret_masking_ggsw_lib");
	}

	add_biv_poly(module, params, acc, acc, phase);

	// The pointer to the last row of the ciphertext vector (B)
	PolyBiv* b_0 = glwe_extract_start_poly(glwe, k);

	// Normalize acc (B) and put it in the result variable
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, b_0, l, (k + 1) * nn, acc, l, nn),
	           "vec_znx_normalize_base2k_p failed in glwe_secret_masking_ggsw_lib");

	status = 0;

cleanup:
	free(as_j);
	free(as_j_dft);
	free(acc);

	return status;
}

int glwe_secret_encrypt_rnx(const MODULE* module, GLWECiphertext* result, const GLWESecretKeyDFT* sk_dft,
                            const PolyUnivRnX* m_univ_rnx)
{
	int status              = -1;
	PolyBiv* biv_phase      = new_biv_poly(result->params);
	PolyUnivRnX* univ_phase = new_univ_rnx(result->params);

	CHECK_CALL(add_normal_random_vec(univ_phase, result->params->nn, m_univ_rnx, 0.0, result->params->sigma),
	           "failed to add the error in glwe encryption");
	CHECK_CALL(univ_rnx_to_biv(result->params, biv_phase, univ_phase),
	           "failed univ to biv conversion in glwe encryption");
	CHECK_CALL(glwe_secret_encrypt_phase(module, result, sk_dft, biv_phase), "masking failed in glwe encryption");

	status = 0;
cleanup:
	free(biv_phase);
	delete_univ_rnx(univ_phase);
	return status;
}

int glwe_secret_encrypt_tnx(const MODULE* module, GLWECiphertext* result, const GLWESecretKeyDFT* sk_dft,
                            const PolyUnivTnX* m_univ_tnx)
{
	int status = -1;

	uint64_t nn = result->params->nn;

	PolyBiv* biv_phase      = new_biv_poly(result->params);
	PolyUnivRnX* err        = new_univ_rnx(result->params);
	PolyUnivTnX* univ_phase = new_univ_tnx(result->params);

	CHECK_CALL(normal_random_vec(err, nn, 0.0, result->params->sigma), "failed to add the error in glwe encryption");
	CHECK_CALL(univ_rnx_to_tnx(result->params, univ_phase, err), "Error in rnx to tnx error conversion for encryption");
	pvda_vec_znx_add(module, univ_phase, 1, nn, univ_phase, 1, nn, m_univ_tnx, 1, nn);

	CHECK_CALL(univ_tnx_to_biv(result->params, biv_phase, univ_phase),
	           "failed univ to biv conversion in glwe encryption");
	CHECK_CALL(glwe_secret_encrypt_phase(module, result, sk_dft, biv_phase), "masking failed in glwe encryption");

	status = 0;
cleanup:
	free(biv_phase);
	delete_univ_rnx(err);
	delete_univ_tnx(univ_phase);
	return status;
}

int glwe_secret_decrypt(const MODULE* module, PolyBiv* res, const GLWESecretKeyDFT* sk_dft, const GLWECiphertext* glwe)
{
	const GLWEParams* params = glwe->params;

	int status = -1;

	// GLWE parameters
	uint64_t nn = params->nn;
	uint64_t k  = params->k;
	uint64_t l  = glwe_params_l(params);

	// Variables
	PolyBiv* acc         = new_biv_poly(params);      // -Sum_j{0,k-1}[sk_j * a_j]
	PolyBivDFT* as_j_dft = new_biv_poly_dft(params);  // DFT(sk_j * a_j)
	PolyBiv* as_j        = new_biv_poly(params);      // sk_j * a_j

	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_demasking_ggsw_lib");
	CHECK_ALLOC(as_j_dft, "as_j_dft's calloc failed in glwe_secret_demasking_ggsw_lib");
	CHECK_ALLOC(as_j, "as_j's calloc failed in glwe_secret_demasking_ggsw_lib");

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++)
	{
		// The j-th component of the secret key in DFT form and the bivGLWE/GLW ciphertext respectively
		PolyUnivDFT* sk_j_univ_dft = glwe_sk_extract_poly_dft(sk_dft, j);
		const PolyUniv* a_j        = glwe_extract_start_poly(glwe, j);

		// Computes DFT(sk_j * a_j)
		pvda_svp_apply_dft(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * nn);

		// Computes sk_j * a_j by inverting the DFT
		CHECK_CALL(pvda_vec_znx_idft(module, as_j, l, as_j_dft, l),
		           "vec_znx_idft_p failed in glwe_secret_demasking_ggsw_lib");

		// And subs it to acc
		pvda_vec_znx_sub(module, acc, l, nn, acc, l, nn, as_j, l, nn);
	}

	// Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
	const PolyBiv* b = glwe_extract_start_poly(glwe, k);

	// acc += b <=> acc = b - sum(sk_j*a_j)
	pvda_vec_znx_add(module, acc, l, nn, acc, l, nn, b, l, (k + 1) * nn);
	//normalize acc into result
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, params->kappa, res, l, nn, acc, l, nn),
	           "vec_znx_normalize_base2k_p failed in glwe_secret_demasking_ggsw_lib");

	status = 0;

cleanup:
	free(as_j);
	free(as_j_dft);
	free(acc);

	return status;
}
