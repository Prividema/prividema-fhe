#include "glwe_ciphertext.h"

#include <string.h>

#include "logger.h"
#include "utils.h"

//! bivGLWE PART (begin)

uint64_t glwe_coef_number(const GLWEParams* params) { return glwe_size(params) * params->nn; }

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

int normalize_glwe(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t nn    = result->params->nn;
	uint64_t k     = result->params->k;
	uint64_t kappa = result->params->kappa;
	uint64_t l     = poly_biv_size(result->params);

	for (uint64_t j = 0; j < k + 1; j++)
		CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, result->vec + j * nn, l, (k + 1) * nn,
		                                         glwe->vec + j * nn, l, (k + 1) * nn),
		           "vec_znx_normalize_base2k_p failed in normalize_glwe");

	status = 0;

cleanup:

	return status;
}

void add_glwe(GLWECiphertext* result, const GLWECiphertext* glwe_lhs, const GLWECiphertext* glwe_rhs)
{
	for (uint64_t p = 0; p < glwe_coef_number(result->params); p++)
		result->vec[p] = glwe_lhs->vec[p] + glwe_rhs->vec[p];
}

int const_mult_glwe(const MODULE* module, GLWECiphertext* result, const PolyUnivDFT* u_dft, const GLWECiphertext* glwe,
                    int do_normalization)
{
	int status = -1;

	// Variables
	VecBivDFT* u_glwe_vec_dft = NULL;

	// bivGLWE parameters
	uint64_t nn = result->params->nn;
	uint64_t k  = result->params->k;
	uint64_t l  = poly_biv_size(result->params);

	// The pointer to DFT(u * glwe)
	u_glwe_vec_dft = malloc(glwe_bytes(result->params));
	CHECK_ALLOC(u_glwe_vec_dft, "u_glwe_dft's malloc failed in const_mult_glwe.");

	// Computes DFT(u * glwe)
	pvda_svp_apply_dft(module, u_glwe_vec_dft, glwe_size(result->params), u_dft, glwe->vec, glwe_size(result->params),
	                   nn);

	// Computes it out of the DFT domain
	CHECK_CALL(
	    pvda_vec_znx_idft(module, result->vec, glwe_size(result->params), u_glwe_vec_dft, glwe_size(result->params)),
	    "vec_znx_idft_p failed in const_mult_glwe");

	if (do_normalization)
		for (uint64_t j = 0; j < k + 1; j++)
			CHECK_CALL(pvda_vec_znx_normalize_base2k(module, result->params->kappa, result->vec + j * nn, l,
			                                         (k + 1) * nn, result->vec + j * nn, l, (k + 1) * nn),
			           "vec_znx_normalize_base2k_p failed in const_mult_glwe");

	status = 0;

cleanup:
	free(u_glwe_vec_dft);

	return status;
}

//! bivGLWE IN DFT PART (begin)

uint64_t glwe_coef_number_dft(const GLWEParams* params) { return glwe_size(params) * params->nn / 2; }

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

int normalize_glwe_dft(const MODULE* module, GLWECiphertextDFT* result_dft, const GLWECiphertextDFT* glwe_dft)
{
	int status = -1;

	// Variables
	VecBiv* glwe_vec = NULL;

	// bivGLWE parameters
	const GLWEParams* params_glwe = result_dft->params;
	uint64_t nn                   = params_glwe->nn;
	uint64_t k                    = params_glwe->k;
	uint64_t kappa                = params_glwe->kappa;
	uint64_t l                    = poly_biv_size(params_glwe);

	// Point to the bivGLWE ciphertext out of the DFT domain
	glwe_vec = malloc(glwe_bytes(params_glwe));

	// Computes the bivGLWE ciphertext out of the DFT domain
	CHECK_CALL(pvda_vec_znx_idft(module, glwe_vec, glwe_size(params_glwe), glwe_dft->vec, glwe_size(params_glwe)),
	           "vec_znx_idft_p failed in normalize_glwe_dft");

	// Normalize the k+1 bivGLWE' elements
	for (uint64_t j = 0; j < k + 1; j++)
		CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, glwe_vec + j * nn, l, (k + 1) * nn, glwe_vec + j * nn,
		                                         l, (k + 1) * nn),
		           "vec_znx_normalize_base2k_p failed in normalize_glwe_dft");

	// Computes the bivGLWE ciphertext in the DFT domain
	pvda_vec_znx_dft(module, result_dft->vec, glwe_size(params_glwe), glwe_vec, glwe_size(params_glwe), nn);

	status = 0;

cleanup:
	free(glwe_vec);

	return status;
}

void add_glwe_dft(GLWECiphertextDFT* result_dft, const GLWECiphertextDFT* glwe_lhs_dft,
                  const GLWECiphertextDFT* glwe_rhs_dft)
{
	for (uint64_t t = 0; t < glwe_coef_number(result_dft->params); t++)
		result_dft->vec[t] = glwe_lhs_dft->vec[t] + glwe_rhs_dft->vec[t];
}

int const_mult_glwe_dft(const MODULE* module, GLWECiphertextDFT* result_dft, const PolyUnivDFT* u_dft,
                        const GLWECiphertextDFT* glwe_dft, int do_normalization)
{
	int status = -1;

	// Variables
	VecBivDFT* u_glwe_vec_dft     = NULL;
	VecBiv* glwe_vec              = NULL;
	VecBiv* result_vec_normalized = NULL;

	// bivGLWE set of parameters
	const GLWEParams* params = result_dft->params;

	// bivGLWE parameters
	uint64_t nn = params->nn;
	uint64_t k  = params->k;
	uint64_t l  = poly_biv_size(params);

	// Point to DFT(u * glwe)
	u_glwe_vec_dft = malloc(glwe_bytes(params));
	CHECK_ALLOC(u_glwe_vec_dft, "u_glwe_dft's malloc failed in const_mult_glwe_dft.");

	// Point to the bivGLWE ciphertext out of the DFT domain
	glwe_vec = malloc(glwe_bytes(params));
	CHECK_ALLOC(glwe_vec, "glwe_vec's malloc failed in const_mult_glwe_dft.");

	CHECK_CALL(pvda_vec_znx_idft(module, glwe_vec, glwe_size(params), glwe_dft->vec, glwe_size(params)),
	           "vec_znx_idft_p failed in const_mult_glwe_dft");

	// Computes DFT(u * glwe)
	pvda_svp_apply_dft(module, result_dft->vec, glwe_size(params), u_dft, glwe_vec, glwe_size(params), nn);

	if (do_normalization)
	{
		// Point to the bivGLWE ciphertext out of the DFT domain
		result_vec_normalized = malloc(glwe_bytes(params));
		CHECK_ALLOC(result_vec_normalized, "result_normalized's malloc failed in const_mult_glwe_dft.");

		// Computes the bivGLWE ciphertext out of the DFT domain
		CHECK_CALL(
		    pvda_vec_znx_idft(module, result_vec_normalized, glwe_size(params), result_dft->vec, glwe_size(params)),
		    "vec_znx_idft_p failed in const_mult_glwe_dft");

		// Normalizes each of the k+1 bivariate polynomials
		for (uint64_t j = 0; j < k + 1; j++)
			CHECK_CALL(pvda_vec_znx_normalize_base2k(module, params->kappa, result_vec_normalized + j * nn, l,
			                                         (k + 1) * nn, result_vec_normalized + j * nn, l, (k + 1) * nn),
			           "vec_znx_normalize_base2k_p failed in const_mult_glwe_dft");

		// Computes the bivGLWE ciphertext in the DFT domain
		pvda_vec_znx_dft(module, result_dft->vec, glwe_size(params), result_vec_normalized, glwe_size(params), nn);
	}

	status = 0;

cleanup:
	free(result_vec_normalized);
	free(glwe_vec);
	free(u_glwe_vec_dft);

	return status;
}

//! COMMON PART (begin)

uint64_t glwe_size(const GLWEParams* params) { return params->n_limbs; }

uint64_t glwe_bytes(const GLWEParams* params)
{
	uint64_t nn = params->nn;
	return glwe_size(params) * nn * sizeof(int64_t);
}

void mult_vec_znx_dft(const MODULE* module, double* result_dft, int64_t result_size, const double* c_dft,
                      int64_t c_size, const double* d_dft, int64_t d_size)
{
	uint64_t nn = pvda_module_extract_nn(module);

	// TODO: we SHOULD offload this to spqlios.
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
