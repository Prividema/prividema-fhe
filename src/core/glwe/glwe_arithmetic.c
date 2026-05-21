#include "glwe_arithmetic.h"

#include <string.h>
#include <sys/types.h>

#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "maths_structures.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

int normalize_glwe(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t k     = result->params->k;
	uint64_t kappa = result->params->kappa;

	for (uint64_t j = 0; j <= k; j++)
	{
		PolyBiv aj_biv  = glwe_extract_poly_view(glwe, j);
		PolyBiv res_biv = glwe_extract_poly_view(result, j);
		CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, &res_biv, &aj_biv),
		           "vec_znx_normalize_base2k_p failed in normalize_glwe");
	}
	status = 0;

cleanup:

	return status;
}

void add_glwe(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe_lhs,
              const GLWECiphertext* glwe_rhs)
{
	PolyBiv lhs_flattened = glwe_flattened_biv(glwe_lhs);
	PolyBiv rhs_flattened = glwe_flattened_biv(glwe_rhs);
	PolyBiv res_flattened = glwe_flattened_biv(result);
	pvda_vec_znx_add(module, &res_flattened, &lhs_flattened, &rhs_flattened);
}

void sub_glwe(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe_lhs,
              const GLWECiphertext* glwe_rhs)
{
	PolyBiv lhs_flattened = glwe_flattened_biv(glwe_lhs);
	PolyBiv rhs_flattened = glwe_flattened_biv(glwe_rhs);
	PolyBiv res_flattened = glwe_flattened_biv(result);
	pvda_vec_znx_sub(module, &res_flattened, &lhs_flattened, &rhs_flattened);
}

void negate_glwe(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe)
{
	uint64_t nn            = result->params->nn;
	PolyBiv glwe_flattened = glwe_flattened_biv(glwe);
	PolyBiv res_flattened  = glwe_flattened_biv(result);
	pvda_vec_znx_negate(module, &res_flattened, &glwe_flattened);
}

int const_mult_glwe(const MODULE* module, GLWECiphertext* result, const PolyUnivDFT* u_dft, const GLWECiphertext* glwe)
{
	int status = -1;

	GLWECiphertextDFT* u_glwe_dft = new_glwe_dft(glwe->params);
	uint64_t nn                   = result->params->nn;

	CHECK_ALLOC(u_glwe_dft, "u_glwe_dft's malloc failed in const_mult_glwe.");

	// Computes DFT(u * glwe)
	PolyBiv glwe_flattened = glwe_flattened_biv(glwe);
	pvda_svp_apply_dft(module, u_glwe_dft->vec, glwe_params_n_limbs(result->params), u_dft, &glwe_flattened);

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

void mult_vec_znx_dft(const MODULE* module, double* result_dft, int64_t result_size, const double* c_dft,
                      int64_t c_size, const double* d_dft, int64_t d_size)
{
	uint64_t nn = pvda_module_extract_nn(module);

	if (c_size <= d_size)
	{
		int64_t smin = c_size < result_size ? c_size : result_size;

		for (uint64_t i = 0; i < smin; i++)
		{
			for (uint64_t j = 0; j < nn / 2; j++)
			{
				// i*N + j corresultponds to the j-th coefficient's index of Re[DFT(c_i)] and Re[DFT(d_i)]
				double c_re = c_dft[i * nn + j];
				double d_re = d_dft[i * nn + j];

				// i*N + j + N/2 corresultponds to the j-th coefficient's index of Im[DFT(c_i)] & Im[DFT(d_i)]
				double c_im = c_dft[i * nn + j + nn / 2];
				double d_im = d_dft[i * nn + j + nn / 2];

				result_dft[i * nn + j]          = c_re * d_re - c_im * d_im;
				result_dft[i * nn + j + nn / 2] = c_re * d_im + c_im * d_re;
			}
		}

		// fill up remaining part with 0's
		double* const dresult_dft = (double*)result_dft;
		memset(dresult_dft + smin * nn, 0, (result_size - smin) * nn * sizeof(double));
	}
	else
		mult_vec_znx_dft(module, result_dft, result_size, d_dft, d_size, c_dft, c_size);
}
