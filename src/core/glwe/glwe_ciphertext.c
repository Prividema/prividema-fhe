#include "glwe_ciphertext.h"

#include <string.h>

#include "logger.h"

//! GLWE PART (begin)

uint64_t glwe_coef_number(const GLWECtParams* params) { return glwe_size(params) * params->N; }

GLWECiphertext* new_glwe(const GLWECtParams* params)
{
	GLWECiphertext* ct = malloc(sizeof(GLWECiphertext));
	if (log_is_null(ct, "ct's malloc failed in new_glwe")) return NULL;

	ct->params = params;

	ct->vec = calloc(glwe_coef_number(params), sizeof(int64_t));
	if (log_is_null(ct->vec, "ct->vec's calloc failed in new_glwe") < 0)
    {
        free(ct);
		return NULL;
    }

	return ct;
}

void delete_glwe(GLWECiphertext* ct)
{
	free(ct->vec);
	free(ct);
}

void normalize_glwe(const MODULE* module, GLWECiphertext* res, const GLWECiphertext* ct_glwe)
{
	// GLWE parameters
	uint64_t N     = res->params->N;
	uint64_t k     = res->params->k;
	uint64_t kappa = res->params->kappa;
	uint64_t l     = poly_biv_size(res->params);

	for (uint64_t j = 0; j < k + 1; j++)
		vec_znx_normalize_base2k_p(module, kappa, res->vec + j * N, l, (k + 1) * N, ct_glwe->vec + j * N, l,
		                           (k + 1) * N);
}

void add_glwe(GLWECiphertext* res, const GLWECiphertext* ct1, const GLWECiphertext* ct2)
{
	for (uint64_t p = 0; p <= glwe_coef_number(res->params); p++)
		res->vec[p] = ct1->vec[p] + ct2->vec[p];
}

int const_mult_glwe(const MODULE* module, GLWECiphertext* res, const PolyUnivDFT* u_dft, const GLWECiphertext* ct, int do_normalization)
{
	// GLWE parameters
	uint64_t N = res->params->N;
	uint64_t k = res->params->k;
	uint64_t l = poly_biv_size(res->params);

	// The pointer to DFT(u * ct)
	VecBivDFT* u_ct_dft = malloc(glwe_bytes(res->params));
	if (log_is_null(u_ct_dft, "u_ct_dft's malloc failed in const_mult_glwe.") < 0) return -1;

	// Computes DFT(u * ct)
	svp_apply_dft_p(module, u_ct_dft, glwe_size(res->params), u_dft, ct->vec, glwe_size(res->params), N);

	// Computes it out of DFT space
	if (vec_znx_idft_p(module, res->vec, glwe_size(res->params), u_ct_dft, glwe_size(res->params)) < 0)
    {
        free(u_ct_dft);
        return log_perror("vec_znx_idft_p failed in const_mult_glwe");
    }

	if (do_normalization)
		for (uint64_t j = 0; j < k + 1; j++)
			vec_znx_normalize_base2k_p(module, res->params->kappa, res->vec + j * N, l, (k + 1) * N, res->vec + j * N,
			                           l, (k + 1) * N);

	free(u_ct_dft);

	return 0;
}

//! GLWE IN DFT PART (begin)

uint64_t glwe_coef_number_dft(const GLWECtParams* params) { return glwe_size(params) * params->N / 2; }

GLWECiphertextDFT* new_glwe_dft(const GLWECtParams* params)
{
	// The GLWE ciphertext in DFT space
	GLWECiphertextDFT* ct_dft = malloc(sizeof(GLWECiphertextDFT));
	if (log_is_null(ct_dft, "ct_dft's malloc failed in new_glwe_dft.") < 0) return NULL;

	// The GLWE parameters
	ct_dft->params = params;

	// The GLWE ciphertext's vector in DFT space
	ct_dft->vec = calloc(glwe_coef_number_dft(params), 2 * sizeof(double));
	if (log_is_null(ct_dft->vec, "ct_dft's calloc failed in new_glwe_dft.") < 0) {
		free(ct_dft);
		return NULL;
	}

	return ct_dft;
}

void delete_glwe_dft(GLWECiphertextDFT* ct)
{
	free(ct->vec);
	free(ct);
}

void add_glwe_dft(GLWECiphertextDFT* res_dft, const GLWECiphertextDFT* ct1_dft, const GLWECiphertextDFT* ct2_dft)
{
    for (uint64_t p = 0; p < glwe_coef_number(res_dft->params); p++)
        res_dft->vec[p] = ct1_dft->vec[p] + ct2_dft->vec[p];
}

int const_mult_glwe_dft(const MODULE* module, GLWECiphertextDFT* res_dft, const PolyUnivDFT* u_dft, const GLWECiphertextDFT* ct_dft,
                        int do_normalization)
{
	// GLWE parameters
	GLWECtParams* params = res_dft->params;
	uint64_t N           = res_dft->params->N;
	uint64_t k           = res_dft->params->k;
	uint64_t l           = poly_biv_size(params);

	// The pointer to DFT(u * ct)
	VecBivDFT* u_ct_dft = malloc(glwe_bytes(params));
	if (log_is_null(u_ct_dft, "u_ct_dft's malloc failed in const_mult_glwe_dft.") < 0) return -1;

	// Computes the GLWE ciphertext out of DFT space
	VecBiv* ct_vec = malloc(glwe_bytes(params));
	if (log_is_null(ct_vec, "ct_vec's malloc failed in const_mult_glwe_dft.") < 0) {
		free(u_ct_dft);
		return -1;
	}
	if (vec_znx_idft_p(module, ct_vec, glwe_size(params), ct_dft->vec, glwe_size(params)) < 0)
    {
        free(u_ct_dft);
        free(ct_vec);

        return log_perror("vec_znx_idft_p failed in const_mult_glwe_dft");
    }

	// Computes DFT(u * ct)
	svp_apply_dft_p(module, res_dft->vec, glwe_size(params), u_dft, ct_vec, glwe_size(params), N);

	free(u_ct_dft);
	free(ct_vec);

	if (do_normalization) {
		// Computes the GLWE ciphertext out of DFT space to normalize it
		VecBiv* res_vec_normalized = malloc(glwe_bytes(params));
		if (log_is_null(res_vec_normalized, "res_vec_normalized's malloc failed in const_mult_glwe_dft.") < 0)
			return -1;
		if (vec_znx_idft_p(module, res_vec_normalized, glwe_size(params), res_dft->vec, glwe_size(params)) < 0)
        {
            free(res_vec_normalized);
            return log_perror("vec_znx_idft_p failed in const_mult_glwe_dft");
        }

		// Normalizes each of the k+1 bivariate polynomials
		for (uint64_t j = 0; j < k + 1; j++)
		{
            if (vec_znx_normalize_base2k_p(module, params->kappa, res_vec_normalized + j * N, l, (k + 1) * N,
			                           res_vec_normalized + j * N, l, (k + 1) * N) < 0)
            {
                free(res_vec_normalized);
                return log_perror("vec_znx_normalize_base2k_p failed in const_mult_glwe_dft");
            }
        }

		// Computes the GLWE ciphertext in DFT space
		vec_znx_dft_p(module, res_dft->vec, glwe_size(params), res_vec_normalized, glwe_size(params), N);

		free(res_vec_normalized);
	}

	return 0;
}

//! COMMON PART (begin)

uint64_t glwe_size(const GLWECtParams* params) { return params->n_limbs; }

uint64_t glwe_bytes(const GLWECtParams* params)
{
	uint64_t N = params->N;
	return glwe_size(params) * N * sizeof(int64_t);
}

void mult_vec_znx_dft(const MODULE* module, double* res_dft, int64_t res_size, const double* c_dft, int64_t c_size,
                      const double* d_dft, int64_t d_size)
{
	uint64_t N = module->nn;

	if (c_size <= d_size) {
		int64_t smin = c_size < res_size ? c_size : res_size;

		for (uint64_t i = 0; i < smin; i++) {
			for (uint64_t j = 0; j < N / 2; j++) {
				// i*N + j corresponds to the j-th coefficient's index of Re[DFT(c_i)] and Re[DFT(d_i)]
				double c_re = c_dft[i * N + j];
				double d_re = d_dft[i * N + j];

				// i*N + j + N/2 corresponds to the j-th coefficient's index of Im[DFT(c_i)] & Im[DFT(d_i)]
				double c_im = c_dft[i * N + j + N / 2];
				double d_im = d_dft[i * N + j + N / 2];

				res_dft[i * N + j]         = c_re * d_re - c_im * d_im;
				res_dft[i * N + j + N / 2] = c_re * d_im + c_im * d_re;
			}
		}

		// fill up remaining part with 0's
		double* const dres_dft = (double*)res_dft;
		memset(dres_dft + smin * N, 0, (res_size - smin) * N * sizeof(double));
	} else
		mult_vec_znx_dft(module, res_dft, res_size, d_dft, d_size, c_dft, c_size);
}