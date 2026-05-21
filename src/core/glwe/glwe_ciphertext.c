#include "glwe_ciphertext.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "maths_structures.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

//! bivGLWE PART (begin)

GLWECiphertext* new_glwe(const GLWEParams* params)
{
	// The bivGLWE does not own the GLWEParams
	GLWECiphertext* glwe = malloc(sizeof(GLWECiphertext));
	CHECK_ALLOC(glwe, "glwe's malloc failed in new_glwe");

	glwe->params = params;

	// Initialize the bivGLWE ciphertext with 0s'
	glwe->vec = calloc(glwe_coef_number(params), sizeof(int64_t));
	CHECK_ALLOC(glwe->vec, "glwe->vec's calloc failed in new_glwe");

	return glwe;
cleanup:

	free(glwe);
	return NULL;
}

void delete_glwe(GLWECiphertext* glwe)
{
	if (!glwe) return;
	free(glwe->vec);
	free(glwe);
}

void glwe_copy(GLWECiphertext* dst, const GLWECiphertext* src)
{
	dst->params = src->params;
	memcpy(dst->vec, src->vec, glwe_coef_number(src->params) * sizeof(int64_t));
}
PolyBiv glwe_extract_poly_view(const GLWECiphertext* glwe_ct, uint64_t pos)
{
	uint64_t nn = glwe_ct->params->nn;
	uint64_t k  = glwe_ct->params->k;
	uint64_t l  = pos == k ? glwe_params_l_b(glwe_ct->params) : glwe_params_l_a(glwe_ct->params);
	PolyBiv ret = {nn, l, (int64_t)((k + 1) * nn), glwe_ct->vec + pos * nn};
	return ret;
}
PolyBiv glwe_flattened_biv(const GLWECiphertext* glwe_ct)
{
	PolyBiv glwe_flattened = {glwe_ct->params->nn, glwe_params_n_limbs(glwe_ct->params), (int64_t)glwe_ct->params->nn,
	                          glwe_ct->vec};
	return glwe_flattened;
}

PolyBivDFT* glwe_extract_start_poly_dft(const GLWECiphertextDFT* glwe_dft, uint64_t pos)
{
	assert(pos >= 0 && pos <= glwe_dft->params->k);
	return glwe_dft->vec + glwe_dft->params->nn * pos;
}

// bivGLWE IN DFT PART (begin)

GLWECiphertextDFT* new_glwe_dft(const GLWEParams* params)
{
	// The bivGLWE ciphertext in the DFT domain
	GLWECiphertextDFT* glwe_dft = malloc(sizeof(GLWECiphertextDFT));
	CHECK_ALLOC(glwe_dft, "glwe_dft's malloc failed in new_glwe_dft.");

	// The bivGLWE parameters
	glwe_dft->params = params;

	// initialize the bivGLWE ciphertext with 0s'
	glwe_dft->vec = calloc(glwe_coef_number_dft(params), 2 * sizeof(double));
	CHECK_ALLOC(glwe_dft->vec, "glwe_dft's calloc failed in new_glwe_dft.");

	return glwe_dft;

cleanup:
	free(glwe_dft);
	return NULL;
}

void delete_glwe_dft(GLWECiphertextDFT* glwe)
{
	if (!glwe) return;
	free(glwe->vec);
	free(glwe);
}

int const_mult_glwe_dft(const MODULE* module, GLWECiphertextDFT* result_dft, const PolyUnivDFT* u_dft,
                        const GLWECiphertextDFT* glwe_dft)
{
	int status               = -1;
	const GLWEParams* params = result_dft->params;
	uint64_t nn              = params->nn;

	VecBiv* glwe_vec       = malloc(glwe_params_bytes(params));
	PolyBiv glwe_flattened = {nn, glwe_params_n_limbs(params), nn, glwe_vec};

	CHECK_ALLOC(glwe_vec, "glwe_vec's malloc failed in const_mult_glwe_dft.");

	CHECK_CALL(pvda_vec_znx_idft(module, &glwe_flattened, glwe_dft->vec, glwe_params_n_limbs(params)),
	           "vec_znx_idft_p failed in const_mult_glwe_dft");

	// Computes DFT(u * glwe)
	//
	pvda_svp_apply_dft(module, result_dft->vec, glwe_params_n_limbs(params), u_dft, &glwe_flattened);

	status = 0;

cleanup:
	free(glwe_vec);

	return status;
}

int glwe_coef_to_dft(const MODULE* module, GLWECiphertextDFT* res_dft, const GLWECiphertext* glwe_ct)
{
	int status = -1;

	PolyBiv glwe_flattened = {glwe_ct->params->nn, glwe_params_n_limbs(glwe_ct->params), glwe_ct->params->nn,
	                          glwe_ct->vec};
	pvda_vec_znx_dft(module, res_dft->vec, glwe_params_n_limbs(glwe_ct->params), &glwe_flattened);

	status = 0;
cleanup:
	return status;
}

int glwe_dft_to_coef(const MODULE* module, GLWECiphertext* res_ct, const GLWECiphertextDFT* glwe_dft)
{
	int status             = -1;
	PolyBiv glwe_flattened = glwe_flattened_biv(res_ct);

	CHECK_CALL(pvda_vec_znx_idft(module, &glwe_flattened, glwe_dft->vec, glwe_params_n_limbs(glwe_dft->params)),
	           "iDFT failed for a GLWE");

	status = 0;
cleanup:
	return status;
}

int glwe_secret_encrypt_phase(const MODULE* module, GLWECiphertext* glwe, const GLWESecretKeyPrepared* sk_prep,
                              const PolyBiv* phase)
{
	int status = -1;

	const GLWEParams* params = (const GLWEParams*)glwe->params;
	// GLWE parameters
	uint64_t nn    = params->nn;
	uint64_t k     = params->k;
	uint64_t kappa = params->kappa;
	uint64_t l_a   = glwe_params_l_a(params);
	uint64_t l_b   = glwe_params_l_b(params);

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
	CHECK_CALL(uniform_random_vec(k * nn, glwe->vec, l_a, (k + 1) * nn, kappa),
	           "A generation failed in glwe_secret_masking_dft");

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	if (k > 1)
	{
		as_j = new_biv_poly(params);  // sk_j * a_j
		CHECK_ALLOC(as_j, "as_j's calloc failed in glwe_secret_masking");

		for (uint64_t j = 0; j < k; j++)
		{
			// The j-th component of the DFT encoding of the secret key
			PolyUnivDFT* sk_j_univ_dft = glwe_prepared_sk_extract_poly_dft(sk_prep, j);

			// Computes DFT(sk_j) * DFT(a_j)
			PolyBiv a_j = glwe_extract_poly_view(glwe, j);
			pvda_svp_apply_dft(module, as_j_dft, l_a, sk_j_univ_dft, &a_j);

			// Undo DFT to retreive sk_j * a_j
			CHECK_CALL(pvda_vec_znx_idft(module, as_j, as_j_dft, l_a),
			           "vec_znx_idft_p failed in glwe_secret_masking_ggsw_lib");

			add_biv_poly(module, params, acc, acc, as_j);
		}
	}
	else
	{
		PolyUnivDFT* sk_j_univ_dft = glwe_prepared_sk_extract_poly_dft(sk_prep, 0);

		PolyBiv a_0 = glwe_extract_poly_view(glwe, 0);
		pvda_svp_apply_dft(module, as_j_dft, l_a, sk_j_univ_dft, &a_0);

		CHECK_CALL(pvda_vec_znx_idft(module, acc, as_j_dft, l_a),
		           "vec_znx_idft_p failed in glwe_secret_masking_ggsw_lib");
	}

	add_biv_poly(module, params, acc, acc, phase);

	// The pointer to the last row of the ciphertext vector (B)
	PolyBiv b = glwe_extract_poly_view(glwe, k);

	// Normalize acc (B) and put it in the result variable
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, &b, acc),
	           "vec_znx_normalize_base2k_p failed in glwe_secret_masking_ggsw_lib");

	status = 0;

cleanup:
	delete_biv(as_j);
	free(as_j_dft);
	delete_biv(acc);

	return status;
}

int glwe_secret_encrypt_rnx(const MODULE* module, GLWECiphertext* result, const GLWESecretKeyPrepared* sk_prep,
                            const PolyUnivRnX* m_univ_rnx)
{
	int status         = -1;
	PolyBiv* biv_phase = new_biv_poly(result->params);
	CHECK_ALLOC(biv_phase, "Bivariate phase allocation failed in GLWE encryption");

	CHECK_CALL(univ_rnx_to_biv(result->params, biv_phase, m_univ_rnx, 0),
	           "failed univ to biv conversion in glwe encryption");
	CHECK_CALL(add_biv_noise(module, result->params, biv_phase, biv_phase), "Noise addition failed in GLWE encryption");

	CHECK_CALL(glwe_secret_encrypt_phase(module, result, sk_prep, biv_phase), "masking failed in glwe encryption");

	status = 0;
cleanup:
	delete_biv(biv_phase);
	return status;
}

int glwe_secret_encrypt_tnx(const MODULE* module, GLWECiphertext* result, const GLWESecretKeyPrepared* sk_prep,
                            const PolyUnivTnX* m_univ_tnx)
{
	int status = -1;

	uint64_t nn        = result->params->nn;
	PolyBiv* biv_phase = new_biv_poly(result->params);
	CHECK_ALLOC(biv_phase, "Bivariate phase allocation failed in GLWE encryption");

	CHECK_CALL(univ_tnx_to_biv(result->params, biv_phase, m_univ_tnx, 0),
	           "failed univ to biv conversion in glwe encryption");
	CHECK_CALL(add_biv_noise(module, result->params, biv_phase, biv_phase), "Noise addition failed in GLWE encryption");
	CHECK_CALL(glwe_secret_encrypt_phase(module, result, sk_prep, biv_phase), "masking failed in glwe encryption");

	status = 0;
cleanup:
	delete_biv(biv_phase);
	return status;
}

int glwe_secret_decrypt(const MODULE* module, PolyBiv* res, const GLWESecretKeyPrepared* sk_prep,
                        const GLWECiphertext* glwe)
{
	const GLWEParams* params = glwe->params;

	int status = -1;

	// GLWE parameters
	uint64_t nn  = params->nn;
	uint64_t k   = params->k;
	uint64_t l_a = glwe_params_l_a(params);
	uint64_t l_b = glwe_params_l_b(params);

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
		PolyUnivDFT* sk_j_univ_dft = glwe_prepared_sk_extract_poly_dft(sk_prep, j);
		//const PolyUniv* a_j        = glwe_extract_start_poly(glwe, j);

		PolyBiv a_j = glwe_extract_poly_view(glwe, j);
		// Computes DFT(sk_j * a_j)

		pvda_svp_apply_dft(module, as_j_dft, l_a, sk_j_univ_dft, &a_j);

		// Computes sk_j * a_j by inverting the DFT
		CHECK_CALL(pvda_vec_znx_idft(module, as_j, as_j_dft, l_a),
		           "vec_znx_idft_p failed in glwe_secret_demasking_ggsw_lib");

		// And subs it to acc
		pvda_vec_znx_sub(module, acc, acc, as_j);
	}

	// Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
	PolyBiv b = glwe_extract_poly_view(glwe, k);

	// acc += b <=> acc = b - sum(sk_j*a_j)
	pvda_vec_znx_add(module, acc, acc, &b);
	//normalize acc into result
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, params->kappa, res, acc),
	           "vec_znx_normalize_base2k_p failed in glwe_secret_demasking_ggsw_lib");

	status = 0;

cleanup:
	delete_biv(as_j);
	free(as_j_dft);
	delete_biv(acc);

	return status;
}
