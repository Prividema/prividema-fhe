#include "bivariate_polynomial.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "logger.h"
#include "rng.h"

//! BIV POLY PART (begin)

uint64_t poly_biv_coef_number(const GLWECtParams* params)
{
	uint64_t N = params->N;
	return poly_biv_size(params) * N;
}

PolyBiv* new_normal_random_biv_poly(const MODULE* module, const GLWECtParams* params)
{
	// GLWE parameters
	uint64_t N     = params->N;
	uint64_t kappa = params->kappa;
	uint64_t l     = params->n_limbs / (params->k + 1);

	// Draws a random univariate polynomial P(X) in Rn[X]
	double* rd_pol_univ = malloc(poly_univ_bytes(params));
	if (log_is_null(rd_pol_univ, "rd_pol_univ's malloc failed.") < 0) return NULL;

	for (int64_t p = 0; p < N; p++)
		if (rand_normal(rd_pol_univ + p, 0.0, params->sigma) < 0) {
			log_perror("rand_normal failed in new_normal_random_biv_poly.");
			free(rd_pol_univ);
			return NULL;
		}

	// Stores the base-2kappa normalized bivariate form Pbiv(X,Y) of P(X)
	PolyBiv* rd_pol = malloc(poly_biv_bytes(params));
	if (log_is_null(rd_pol, "rd_pol's malloc failed in new_normal_random_biv_poly.") < 0) {
		free(rd_pol_univ);
		return NULL;
	}

	// For each, (p,i) in [0,N-1]x[0,l-1], Pbiv_p_i = centered(floor(P_p * (2^kappa)^i))
	// Where centered(_) is in [-2^(kappa-1) ; 2^(kappa-1) - 1]
	int64_t mask = (1LL << kappa) - 1;
	for (int64_t p = 0; p < N; p++) {
		// For each p, we substract (2^kappa)/2 * (2^kappa)^(-i) to P_p
		for (int64_t i = 1; i < l + 1; i++) {
			rd_pol_univ[p] += ldexp(1.0, kappa - 1 - kappa * i);
		}

		if (rd_pol_univ[p] >= 0) {
			for (int64_t i = 1; i < l + 1; i++) {
				rd_pol[(i - 1) * N + p] = (((int64_t)ldexp(rd_pol_univ[p], i * kappa)) & mask) - (1LL << (kappa - 1));
			}
		} else {
			rd_pol_univ[p] -= floor(rd_pol_univ[p]);
			for (int64_t i = 1; i < l + 1; i++) {
				rd_pol[(i - 1) * N + p] = (((int64_t)ldexp(rd_pol_univ[p], i * kappa)) & mask) - (1LL << (kappa - 1));
			}
		}
	}

	free(rd_pol_univ);

	return rd_pol;
}

PolyBiv* new_uniform_random_biv_poly(const MODULE* module, const GLWECtParams* params, int64_t precision)
{
	PolyBiv* pol = malloc(poly_biv_bytes(params));
	if (log_is_null(pol, "pol's malloc failed in new_uniform_random_biv_poly.") < 0) return NULL;

	for (int64_t i = 1; i < precision + 1; i++)
		for (int64_t p = 0; p < params->N; p++)
			if (rand_uniform(pol + (i - 1) * params->N + p, params->kappa) < 0) {
				log_perror("rand_uniform failed in new_uniform_random_biv_poly.");
				free(pol);
				return NULL;
			}

	return pol;
}

void add_biv_poly(const GLWECtParams* params, PolyBiv* res, int64_t res_sl, const PolyBiv* a, int64_t a_sl, const PolyBiv* b,
                  int64_t b_sl)
{
	for (int64_t i = 1; i <= poly_biv_size(params); i++)
		for (int64_t p = 0; p < params->N; p++)
			res[p + (i - 1) * res_sl] = a[p + (i - 1) * a_sl] + b[p + (i - 1) * b_sl];
}

//! BIV POLY IN DFT PART (begin)

uint64_t poly_biv_coef_number_dft(const GLWECtParams* params)
{
	uint64_t N = params->N;
	return (poly_biv_size(params) * N) / 2;
}

PolyBivDFT* new_normal_random_biv_poly_dft(const MODULE* module, const GLWECtParams* params)
{
	// Base-2Kappa normalized bivariate polynomial
	PolyBiv* rd_pol = new_normal_random_biv_poly(module, params);
	if (log_is_null(rd_pol, "new_normal_random_biv_poly failed in new_normal_biv_poly_dft.") < 0) return NULL;

	// Base-2Kappa normalized bivariate polynomial in DFT space
	PolyBivDFT* rd_pol_dft = malloc(poly_biv_bytes(params));
	if (log_is_null(rd_pol_dft, "rd_pol_dft's malloc failed in new_normal_random_biv_poly_dft.") < 0) {
		free(rd_pol);
		return NULL;
	}

	// Then compute in DFT space
	vec_znx_dft_p(module, rd_pol_dft, poly_biv_size(params), rd_pol, poly_biv_size(params), params->N);

	free(rd_pol);

	return rd_pol_dft;
}

PolyBivDFT* new_uniform_random_biv_poly_dft(const MODULE* module, const GLWECtParams* params, int64_t precision)
{
	PolyBiv* pol = new_uniform_random_biv_poly(module, params, precision);
	if (log_is_null(pol, "pol's malloc failed in new_uniform_random_biv_poly_dft.") < 0) return NULL;

	PolyBivDFT* pol_dft = malloc(poly_biv_bytes(params));
	if (log_is_null(pol_dft, "pol_dft's malloc failed in new_uniform_random_biv_poly_dft.") < 0) {
		free(pol);
		return NULL;
	}

	vec_znx_dft_p(module, pol_dft, poly_biv_size(params), pol, poly_biv_size(params), params->N);

	free(pol);

	return pol_dft;
}

void add_biv_poly_dft(const GLWECtParams* params, PolyBivDFT* res, int64_t res_sl, const PolyBivDFT* a, int64_t a_sl, const PolyBivDFT* b,
                      int64_t b_sl)
{
	for (int64_t i = 1; i <= poly_biv_size(params); i++)
		for (int64_t p = 0; p < params->N; p++)
			res[p + (i - 1) * res_sl] = a[p + (i - 1) * a_sl] + b[p + (i - 1) * b_sl];
}

//! COMMON PART (begin)

uint64_t poly_biv_bytes(const GLWECtParams* params) { return poly_biv_coef_number(params) * sizeof(int64_t); }

uint64_t poly_biv_size(const GLWECtParams* params) { return params->n_limbs / (params->k + 1); }

uint64_t poly_univ_bytes(const GLWECtParams* params)
{
	uint64_t N = params->N;
	return N * sizeof(int64_t);
}

void biv_to_univ(const GLWECtParams* params, double* res_univ, const PolyBiv* pol_biv)
{
	// GLWE parameters
	uint64_t N     = params->N;
	uint64_t kappa = params->kappa;
	uint64_t l     = poly_biv_size(params);

	// res_univ(X^p) = Sum_i{1,l}[poly(X^p, Y^i) * 2^(-kappa*i)]
	for (int64_t i = 1; i <= l; i++)
		for (int64_t p = 0; p < N; p++) res_univ[p] += ldexp((double)pol_biv[(i - 1) * N + p], -i * kappa);
}

int univ_to_biv(const GLWECtParams* params, PolyBiv* res, const double* pol_univ)
{
	// GLWE parameters
	uint64_t N     = params->N;
	uint64_t kappa = params->kappa;
	uint64_t l     = poly_biv_size(params);

	// Fills each pol_biv(X^p, Y^i) with the pol_univ's decomposition coefficients of  in [-2^(kappa* - 1) ; 2^(kappa -
	// 1) - 1]
	int64_t mask         = (1LL << kappa) - 1;

	double* tmp_pol_univ = malloc(poly_biv_bytes(params));
	if (log_is_null(tmp_pol_univ, "tmp_pol_univ's malloc failed in univ_to_biv.") < 0) return -1;

	for (int64_t p = 0; p < N; p++) {
		// For each p, we substract (2^kappa)/2 * (2^kappa)^(-i) to pol_univ[p]
		tmp_pol_univ[p] = pol_univ[p];
		for (int64_t i = 1; i <= l; i++) {
			tmp_pol_univ[p] += ldexp(1.0, kappa - 1 - kappa * i);
		}

		if (tmp_pol_univ[p] >= 0) {
			for (int64_t i = 1; i <= l; i++) {
				res[(i - 1) * N + p] = (((int64_t)ldexp(tmp_pol_univ[p], i * kappa)) & mask) - (1LL << (kappa - 1));
			}
		} else {
			tmp_pol_univ[p] -= floor(tmp_pol_univ[p]);
			for (int64_t i = 1; i <= l; i++) {
				res[(i - 1) * N + p] = (((int64_t)ldexp(tmp_pol_univ[p], i * kappa)) & mask) - (1LL << (kappa - 1));
			}
		}
	}

	free(tmp_pol_univ);

	return 0;
}