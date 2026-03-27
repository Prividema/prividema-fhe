#include "glwe.h"

#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "logger.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

//! bivGLWE PART (begin)
int add_mult(const MODULE* module, const GLWEParams* params, PolyBiv* res, const VecBiv* glwe,
             const GLWESecretKeyDFT* sk_dft)
{
	int status = -1;
	// GLWE parameters
	uint64_t nn = params->nn;
	uint64_t k  = params->k;
	uint64_t l  = glwe_params_l(params);

	PolyBivDFT* as_j_dft = new_biv_poly_dft(params);  //DFT(sk_j * a_j)
	PolyBiv* as_j        = new_biv_poly(params);      // sk_j * a_j

	CHECK_ALLOC(as_j_dft, "as_j_dft's malloc failed in add_mult");
	CHECK_ALLOC(as_j, "as_j's malloc failed in add_mult");

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++)
	{
		// The j-th component of resp. the secret key and the bivGLWE ciphertext
		PolyUnivDFT* sk_j_univ_dft = glwe_sk_extract_poly_dft(sk_dft, j);

		const PolyBiv* a_j = glwe + j * nn;

		// Computes DFT(sk_j * a_j)
		pvda_svp_apply_dft(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * nn);

		// Invert DFT to get sk_j * a_j
		CHECK_CALL(pvda_vec_znx_idft(module, as_j, l, as_j_dft, l), "vec_znx_idft_p failed in add_mult");

		// Computes acc = acc - sk_j * a_j
		for (uint64_t p = 0; p < nn * l; p++) res[p] += as_j[p];
	}
	status = 0;
cleanup:

	free(as_j_dft);
	free(as_j);

	return status;
}

int sub_mult(const MODULE* module, const GLWEParams* params, PolyBiv* res, const VecBiv* ct,
             const GLWESecretKeyDFT* sk_dft)
{
	int status = -1;
	// GLWE parameters
	uint64_t nn = params->nn;
	uint64_t k  = params->k;
	uint64_t l  = glwe_params_l(params);

	PolyBivDFT* as_j_dft = new_biv_poly_dft(params);  //DFT(sk_j * a_j)
	PolyBiv* as_j        = new_biv_poly(params);      // sk_j * a_j

	CHECK_ALLOC(as_j_dft, "as_j_dft's malloc failed in sub_mult");
	CHECK_ALLOC(as_j, "as_j's malloc failed in sub_mult");

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++)
	{
		// The j-th component of resp. the secret key and the bivGLWE ciphertext
		PolyUnivDFT* sk_j_univ_dft = glwe_sk_extract_poly_dft(sk_dft, j);
		const PolyBiv* a_j         = ct + j * nn;

		// Computes DFT(sk_j * a_j)
		pvda_svp_apply_dft(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * nn);

		// Invert DFT to get sk_j * a_j
		CHECK_CALL(pvda_vec_znx_idft(module, as_j, l, as_j_dft, l), "vec_znx_idft_p failed in sub_mult");

		// Computes acc = acc - sk_j * a_j
		for (uint64_t p = 0; p < nn * l; p++)
		{
			res[p] -= as_j[p];
		}
	}
	status = 0;
cleanup:

	free(as_j_dft);
	free(as_j);

	return status;
}

int glwe_secret_masking(const MODULE* module, GLWECiphertext* glwe, const GLWESecretKeyDFT* sk_dft,
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
	PolyBiv* as_j        = new_biv_poly(params);      // sk_j * a_j

	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_masking");
	CHECK_ALLOC(as_j_dft, "as_j_dft's calloc failed in glwe_secret_masking");
	CHECK_ALLOC(as_j, "as_j's calloc failed in glwe_secret_masking");

	// Draws uniformly in Zn[X,Y] the ajs'
	CHECK_CALL(uniform_random_vec(k * nn, glwe->vec, l, (k + 1) * nn, kappa),
	           "A generation failed in glwe_secret_masking_dft");

	// Computes Sum_j{0,k-1}[sk_j * a_j]
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

int glwe_secret_demasking(const MODULE* module, PolyBiv* res, const GLWESecretKeyDFT* sk_dft,
                          const GLWECiphertext* glwe)
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

//! bivGLWE IN DFT SPACE PART (begin)

int glwe_secret_masking_dft(const MODULE* module, GLWECiphertextDFT* glwe_dft, const GLWESecretKeyDFT* sk_dft,
                            const PolyBivDFT* phase_dft)
{
	int status = -1;
	// GLWE parameters
	const GLWEParams* params = glwe_dft->params;
	uint64_t k               = params->k;
	uint64_t nn              = params->nn;
	uint64_t kappa           = params->kappa;
	uint64_t l               = glwe_params_l(params);

	GLWECiphertext* glwe_ct = new_glwe(params);
	// acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
	PolyBiv* acc = new_biv_poly(params);
	CHECK_ALLOC(glwe_ct, "calloc failed in glwe_secret_masking_dft");
	CHECK_ALLOC(acc, "calloc failed in glwe_secret_masking_dft");

	CHECK_CALL(uniform_random_vec(k * nn, glwe_ct->vec, l, (k + 1) * nn, kappa),
	           "Random vec generation failed in glwe_secret_masking_dft");

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	CHECK_CALL(add_mult(module, params, acc, glwe_ct->vec, sk_dft), "add_mult failed in glwe_secret_masking_dft.");

	// The pointer to limb_0(b)
	PolyBiv* b_0 = glwe_extract_start_poly(glwe_ct, k);

	// For each i in {0,l} limb_i(b) = limb_i(acc) = Sum_j{0,k-1}[sk_j * limb_i(a_j)]
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, b_0, l, (k + 1) * nn, acc, l, nn),
	           "vec_znx_normalize_base2k_p failed in glwe_secret_masking_dft");

	// Computes the bivGLWE ciphertext in the DFT domain
	glwe_coef_to_dft(module, glwe_dft, glwe_ct);

	// Add the phase to the result ciphertext's b
	// TODO:: vec_add
	for (uint64_t i = 0; i < l; i++)
	{
		for (uint64_t p = 0; p < nn; p++)
		{
			glwe_dft->vec[i * (k + 1) * nn + k * nn + p] += phase_dft[i * nn + p];
		}
	}
	/*
	  pvda_vec_rnx_add(module, glwe_dft->vec + k * nn, l, (k + 1) * nn, glwe_dft->vec + k * nn, l, (k + 1) * nn,
	                   phase_dft, l, nn);
	*/

	status = 0;
cleanup:
	delete_glwe(glwe_ct);
	free(acc);

	return status;
}

int glwe_secret_demasking_dft(const MODULE* module, PolyBiv* res, const GLWESecretKeyDFT* sk_dft,
                              const GLWECiphertextDFT* glwe_dft)
{
	int status = -1;
	// GLWE parameters
	const GLWEParams* params = glwe_dft->params;
	uint64_t nn              = params->nn;
	uint64_t k               = params->k;
	uint64_t l               = glwe_params_l(params);

	// Computes the input ciphertext out of the DFT domain
	GLWECiphertext* glwe_ct = new_glwe(glwe_dft->params);
	PolyBiv* acc            = new_biv_poly(glwe_dft->params);
	CHECK_ALLOC(glwe_ct, "calloc failed in glwe_secret_masking_dft");
	CHECK_ALLOC(acc, "calloc failed in glwe_secret_masking_dft");

	CHECK_CALL(
	    pvda_vec_znx_idft(module, glwe_ct->vec, glwe_total_nlimbs(params), glwe_dft->vec, glwe_total_nlimbs(params)),
	    "idft failed in glwe_secret_demasking_dft");

	CHECK_CALL(sub_mult(module, params, acc, glwe_ct->vec, sk_dft), "sub_mult failed in glwe_secret_masking_dft");

	// Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
	PolyBiv* b = glwe_extract_start_poly(glwe_ct, k);
	pvda_vec_znx_add(module, acc, l, nn, acc, l, nn, b, l, (k + 1) * nn);

	// The phase in Zn[X,Y]
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, params->kappa, res, l, nn, acc, l, nn),
	           "normalization failed in glwe_secret_demasking_dft");

	status = 0;

cleanup:
	free(acc);
	delete_glwe(glwe_ct);

	return status;
}
