#include "glwe_arithmetic.h"

#include <string.h>
#include <sys/types.h>

#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

int normalize_glwe(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t nn    = result->params->nn;
	uint64_t k     = result->params->k;
	uint64_t kappa = result->params->kappa;
	uint64_t l_a   = glwe_params_l_a(result->params);
	uint64_t l_b   = glwe_params_l_b(result->params);

	for (uint64_t j = 0; j < k; j++)
		CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, result->vec + j * nn, l_a, (k + 1) * nn,
		                                         glwe->vec + j * nn, l_a, (k + 1) * nn),
		           "vec_znx_normalize_base2k_p failed in normalize_glwe");
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, result->vec + k * nn, l_b, (k + 1) * nn, glwe->vec + k * nn,
	                                         l_b, (k + 1) * nn),
	           "vec_znx_normalize_base2k_p failed in normalize_glwe");
	status = 0;

cleanup:

	return status;
}

void add_glwe(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe_lhs,
              const GLWECiphertext* glwe_rhs)
{
	uint64_t nn = result->params->nn;
	pvda_vec_znx_add(module, result->vec, glwe_params_n_limbs(result->params), nn, glwe_lhs->vec,
	                 glwe_params_n_limbs(glwe_lhs->params), nn, glwe_rhs->vec, glwe_params_n_limbs(glwe_rhs->params),
	                 nn);
}

void sub_glwe(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe_lhs,
              const GLWECiphertext* glwe_rhs)
{
	uint64_t nn = result->params->nn;
	pvda_vec_znx_sub(module, result->vec, glwe_params_n_limbs(result->params), nn, glwe_lhs->vec,
	                 glwe_params_n_limbs(glwe_lhs->params), nn, glwe_rhs->vec, glwe_params_n_limbs(glwe_rhs->params),
	                 nn);
}

void negate_glwe(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe)
{
	uint64_t nn = result->params->nn;
	pvda_vec_znx_negate(module, result->vec, glwe_params_n_limbs(result->params), nn, glwe->vec,
	                    glwe_params_n_limbs(glwe->params), nn);
}

int const_mult_glwe(const MODULE* module, GLWECiphertext* result, const PolyUnivDFT* u_dft, const GLWECiphertext* glwe)
{
	int status = -1;

	GLWECiphertextDFT* u_glwe_dft = new_glwe_dft(glwe->params);
	uint64_t nn                   = result->params->nn;

	CHECK_ALLOC(u_glwe_dft, "u_glwe_dft's malloc failed in const_mult_glwe.");

	// Computes DFT(u * glwe)
	pvda_svp_apply_dft(module, u_glwe_dft->vec, glwe_params_n_limbs(result->params), u_dft, glwe->vec,
	                   glwe_params_n_limbs(result->params), nn);

	// Computes it out of the DFT domain
	glwe_dft_to_coef(module, result, u_glwe_dft);

	status = 0;
cleanup:
	delete_glwe_dft(u_glwe_dft);

	return status;
}
void add_glwe_dft(GLWECiphertextDFT* result_dft, const GLWECiphertextDFT* glwe_lhs_dft,
                  const GLWECiphertextDFT* glwe_rhs_dft)
{
	//TODO: It would be good to move this to SPQLIOS/other backend for faster addition
	// Not being done right now due to header namespace clash if we blindly include
	// spqlios rnx functions.
	for (uint64_t t = 0; t < glwe_coef_number(result_dft->params); t++)
		result_dft->vec[t] = glwe_lhs_dft->vec[t] + glwe_rhs_dft->vec[t];
}
