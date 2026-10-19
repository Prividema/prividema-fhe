#include "glwe_arithmetic.h"

#include <string.h>
#include <sys/types.h>

#include "backend_arithmetic.h"
#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "maths_structures.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

int glwe_normalize(const PvdaBackend* module, GLWECiphertext* result, const GLWECiphertext* glwe)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t k     = result->params->k;
	uint64_t kappa = result->params->kappa;

	for (uint64_t j = 0; j <= k; j++)
	{
		PolyBiv aj_biv  = glwe_extract_poly_view(glwe, j);
		PolyBiv res_biv = glwe_extract_poly_view(result, j);
		CHECK_CALL(biv_normalize(module, kappa, &res_biv, &aj_biv),
		           "vec_znx_normalize_base2k_p failed in glwe_normalize");
	}
	status = 0;

cleanup:

	return status;
}

int glwe_add(const PvdaBackend* module, GLWECiphertext* result, const GLWECiphertext* glwe_lhs,
              const GLWECiphertext* glwe_rhs)
{
	ZnxVec lhs_flattened = glwe_as_znx_vec(glwe_lhs);
	ZnxVec rhs_flattened = glwe_as_znx_vec(glwe_rhs);
	ZnxVec res_flattened = glwe_as_znx_vec(result);
	return pvda_vec_znx_add(module, &res_flattened, &lhs_flattened, &rhs_flattened);
}

int glwe_sub(const PvdaBackend* module, GLWECiphertext* result, const GLWECiphertext* glwe_lhs,
              const GLWECiphertext* glwe_rhs)
{
	ZnxVec lhs_flattened = glwe_as_znx_vec(glwe_lhs);
	ZnxVec rhs_flattened = glwe_as_znx_vec(glwe_rhs);
	ZnxVec res_flattened = glwe_as_znx_vec(result);
	return pvda_vec_znx_sub(module, &res_flattened, &lhs_flattened, &rhs_flattened);
}

int glwe_negate(const PvdaBackend* module, GLWECiphertext* result, const GLWECiphertext* glwe)
{
	uint64_t nn           = result->params->nn;
	ZnxVec glwe_flattened = glwe_as_znx_vec(glwe);
	ZnxVec res_flattened  = glwe_as_znx_vec(result);
	return pvda_vec_znx_negate(module, &res_flattened, &glwe_flattened);
}

int glwe_mult_const(const PvdaBackend* module, GLWECiphertext* result, const BE_PolySvpPrep* u_dft,
                    const GLWECiphertext* glwe)
{
	int status = -1;

	GLWECiphertextDFT* u_glwe_dft = new_glwe_dft(module, glwe->params);
	uint64_t nn                   = result->params->nn;

	CHECK_ALLOC(u_glwe_dft, "u_glwe_dft's malloc failed in glwe_mult_const.");

	// Computes DFT(u * glwe)
	ZnxVec glwe_flattened = glwe_as_znx_vec(glwe);
	pvda_svp_apply_dft(module, u_glwe_dft->vec, glwe_params_n_limbs(result->params), u_dft, &glwe_flattened);

	// Computes it out of the DFT domain
	glwe_dft_to_coef(module, result, u_glwe_dft);

	status = 0;
cleanup:
	delete_glwe_dft(u_glwe_dft);

	return status;
}

int glwe_dft_add(const PvdaBackend* module, GLWECiphertextDFT* result_dft, const GLWECiphertextDFT* glwe_lhs_dft,
                  const GLWECiphertextDFT* glwe_rhs_dft)
{
	// Computes result = lhs + rhs in the DFT domain, one polynomial per limb
	return pvda_vec_znx_dft_add(module, result_dft->vec, glwe_params_n_limbs(result_dft->params), glwe_lhs_dft->vec,
								glwe_params_n_limbs(glwe_lhs_dft->params), glwe_rhs_dft->vec,
								glwe_params_n_limbs(glwe_rhs_dft->params));
}


int glwe_dft_mult_const(const PvdaBackend* module, GLWECiphertextDFT* result_dft, const BE_PolySvpPrep* u_dft,
                        const GLWECiphertextDFT* glwe_dft)
{
	// Computes DFT(u * glwe)
	//
	pvda_svp_apply_dft_to_dft(module, result_dft->vec, glwe_params_n_limbs(result_dft->params), u_dft, glwe_dft->vec,
	                          glwe_params_n_limbs(glwe_dft->params));
	return 0;
}

