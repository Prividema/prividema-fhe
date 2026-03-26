#include "bivariate_polynomial.h"

#include <math.h>

#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

//! BIV POLY PART (begin)

uint64_t poly_biv_coef_number(const GLWEParams* params_glwe)
{
	uint64_t nn = params_glwe->nn;
	return poly_biv_size(params_glwe) * nn;
}

PolyBiv* new_biv_poly(const GLWEParams* params_glwe)
{
	PolyBiv* pol = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	CHECK_ALLOC(pol, "pol's malloc failed in new_biv_poly");
	return pol;
cleanup:
	return NULL;
}

int normal_random_biv_poly(const GLWEParams* params_glwe, PolyBiv* result)
{
	int status = -1;

	// Variables
	PolyUnivRnX* rd_pol_univ = new_univ_rnx(params_glwe);
	CHECK_ALLOC(rd_pol_univ, "rd_pol_univ's malloc failed.");
	// bivGLWE parameters
	uint64_t nn    = params_glwe->nn;
	uint64_t kappa = params_glwe->kappa;
	uint64_t l     = params_glwe->n_limbs / (params_glwe->k + 1);

	// Draws a random univariate polynomial P(X) in Rn[X]
	CHECK_CALL(normal_random_vec(rd_pol_univ, params_glwe->nn, 0.0, params_glwe->sigma),
	           "random normal vec generation failed");

	univ_to_biv(params_glwe, result, rd_pol_univ);

	status = 0;

cleanup:
	delete_univ_rnx(rd_pol_univ);

	return status;
}

int uniform_random_biv_poly(const GLWEParams* params_glwe, PolyBiv* result, int64_t precision)
{
	int status = -1;

	for (uint64_t i = 1; i <= precision; i++)
		for (uint64_t p = 0; p < params_glwe->nn; p++)
			CHECK_CALL(rand_uniform(result + (i - 1) * params_glwe->nn + p, params_glwe->kappa),
			           "rand_uniform failed in uniform_random_biv_poly.");

	status = 0;

cleanup:

	return status;
}

void add_biv_poly(const GLWEParams* params_glwe, PolyBiv* res, int64_t res_sl, const PolyBiv* a, int64_t a_sl,
                  const PolyBiv* b, int64_t b_sl)
{
	for (uint64_t i = 0; i < poly_biv_size(params_glwe); i++)
		for (uint64_t p = 0; p < params_glwe->nn; p++) res[p + i * res_sl] = a[p + i * a_sl] + b[p + i * b_sl];
}

//! BIV POLY IN DFT PART (begin)

uint64_t poly_biv_coef_number_dft(const GLWEParams* params_glwe)
{
	uint64_t nn = params_glwe->nn;
	return (poly_biv_size(params_glwe) * nn) / 2;
}

PolyBivDFT* new_biv_poly_dft(const GLWEParams* params_glwe)
{
	PolyBivDFT* pol_dft = calloc(poly_biv_coef_number_dft(params_glwe), 2 * sizeof(double));
	CHECK_ALLOC(pol_dft, "pol_dft's malloc failed in new_biv_poly");
	return pol_dft;
cleanup:
	return NULL;
}

int normal_random_biv_poly_dft(const MODULE* module, const GLWEParams* params_glwe, PolyBivDFT* result_dft)
{
	int status = -1;

	// Variables
	PolyBiv* rd_pol = NULL;

	// Base-2Kappa normalized bivariate polynomial
	rd_pol = new_biv_poly(params_glwe);
	CHECK_ALLOC(rd_pol, "rd_pol malloc failed in normal_random_biv_poly_dft");
	CHECK_CALL(normal_random_biv_poly(params_glwe, rd_pol), "normal_random_biv_poly failed in normal_biv_poly_dft.");

	// Then compute in the DFT domain
	pvda_vec_znx_dft(module, result_dft, poly_biv_size(params_glwe), rd_pol, poly_biv_size(params_glwe),
	                 params_glwe->nn);

	status = 0;

cleanup:
	free(rd_pol);

	return status;
}

int uniform_random_biv_poly_dft(const MODULE* module, const GLWEParams* params_glwe, PolyBivDFT* result_dft,
                                int64_t precision)
{
	int status = -1;

	// Variables
	PolyBiv* pol = NULL;

	// Uniformly drawn bivariate polynomial
	CHECK_CALL(uniform_random_biv_poly(params_glwe, pol, precision),
	           "pol's malloc failed in uniform_random_biv_poly_dft.");

	// Computes bivariate polynomial in the DFT domain
	pvda_vec_znx_dft(module, result_dft, poly_biv_size(params_glwe), pol, poly_biv_size(params_glwe), params_glwe->nn);

	status = 0;

cleanup:
	free(pol);

	return status;
}

void add_biv_poly_dft(const GLWEParams* params_glwe, PolyBivDFT* res, int64_t res_sl, const PolyBivDFT* a, int64_t a_sl,
                      const PolyBivDFT* b, int64_t b_sl)
{
	for (uint64_t i = 0; i < poly_biv_size(params_glwe); i++)
		for (uint64_t p = 0; p < params_glwe->nn; p++) res[p + i * res_sl] = a[p + i * a_sl] + b[p + i * b_sl];
}

//! COMMON PART (begin)

uint64_t poly_biv_bytes(const GLWEParams* params_glwe) { return poly_biv_coef_number(params_glwe) * sizeof(int64_t); }

uint64_t poly_biv_size(const GLWEParams* params_glwe) { return params_glwe->n_limbs / (params_glwe->k + 1); }

void biv_to_univ(const GLWEParams* params_glwe, double* res_univ, const PolyBiv* pol_biv)
{
	// bivGLWE parameters
	uint64_t nn    = params_glwe->nn;
	uint64_t kappa = params_glwe->kappa;
	uint64_t l     = poly_biv_size(params_glwe);

	// res_univ(X^p) = Sum_i{1,l}[poly(X^p, Y^i) * 2^(-kappa*i)]
	// TODO: slow polynomial evaulation,
	for (uint64_t p = 0; p < nn; p++)
	{
		res_univ[p] = 0;
		for (uint64_t i = 1; i <= l; i++) res_univ[p] += ldexp((double)pol_biv[(i - 1) * nn + p], -i * kappa);
	}
}

int univ_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const double* pol_univ)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t nn    = params_glwe->nn;
	uint64_t kappa = params_glwe->kappa;
	uint64_t l     = poly_biv_size(params_glwe);

	// Fills each pol_biv(X^p, Y^i) with the pol_univ's decomposition coefficients of  in [-2^(kappa* - 1) ; 2^(kappa -
	// 1) - 1]
	int64_t mask = (1LL << kappa) - 1;

	for (uint64_t p = 0; p < nn; p++)
	{
		// For each p, we substract (2^kappa)/2 * (2^kappa)^(-i) to pol_univ[p]
		double tmp = pol_univ[p];
		for (uint64_t i = 1; i <= l; i++)
		{
			tmp += ldexp(1.0, kappa - 1 - kappa * i);
		}

		if (tmp < 0) tmp -= floor(tmp);
		for (uint64_t i = 1; i <= l; i++)
		{
			res[(i - 1) * nn + p] = (((int64_t)ldexp(tmp, i * kappa)) & mask) - (1LL << (kappa - 1));
		}
	}

	status = 0;

cleanup:

	return status;
}

int biv_coefs_to_dft(const MODULE* module, const GLWEParams* params_glwe, PolyBivDFT* res_dft, const PolyBiv* a)
{
	uint64_t nn = params_glwe->nn;
	uint64_t l  = poly_biv_size(params_glwe);
	pvda_vec_znx_dft(module, res_dft, l, a, l, nn);
	return 1;
}

int biv_dft_to_coefs(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBivDFT* a_dft)
{
	uint64_t nn = params_glwe->nn;
	uint64_t l  = poly_biv_size(params_glwe);
	return pvda_vec_znx_idft(module, res, l, a_dft, l);
}
