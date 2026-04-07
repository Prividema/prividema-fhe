#include "bivariate_polynomial.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>

#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

//! BIV POLY PART (begin)

uint64_t poly_biv_coef_number(const GLWEParams* params_glwe)
{
	uint64_t nn = params_glwe->nn;
	return glwe_params_l(params_glwe) * nn;
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

	uint64_t nn    = params_glwe->nn;
	uint64_t kappa = params_glwe->kappa;
	uint64_t l     = params_glwe->n_limbs / (params_glwe->k + 1);

	PolyUnivRnX* rd_pol_univ = new_univ_rnx(params_glwe);
	CHECK_ALLOC(rd_pol_univ, "rd_pol_univ's malloc failed.");

	//TODO: generate directly? Use tnx? Seems like it could use improvement
	CHECK_CALL(normal_random_vec(rd_pol_univ, params_glwe->nn, 0.0, params_glwe->sigma),
	           "random normal vec generation failed");

	univ_rnx_to_biv(params_glwe, result, rd_pol_univ);

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

void add_biv_poly(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBiv* a, const PolyBiv* b)
{
	uint64_t nn = params_glwe->nn;
	uint64_t l  = glwe_params_l(params_glwe);
	pvda_vec_znx_add(module, res, l, nn, a, l, nn, b, l, nn);
}

//! BIV POLY IN DFT PART (begin)

uint64_t poly_biv_coef_number_dft(const GLWEParams* params_glwe)
{
	uint64_t nn = params_glwe->nn;
	return (glwe_params_l(params_glwe) * nn) / 2;
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

	PolyBiv* rd_pol = new_biv_poly(params_glwe);
	CHECK_ALLOC(rd_pol, "rd_pol malloc failed in normal_random_biv_poly_dft");

	CHECK_CALL(normal_random_biv_poly(params_glwe, rd_pol), "normal_random_biv_poly failed in normal_biv_poly_dft.");
	biv_coefs_to_dft(module, params_glwe, result_dft, rd_pol);

	status = 0;
cleanup:
	free(rd_pol);

	return status;
}

int uniform_random_biv_poly_dft(const MODULE* module, const GLWEParams* params_glwe, PolyBivDFT* result_dft,
                                int64_t precision)
{
	int status = -1;

	PolyBiv* pol = NULL;

	CHECK_CALL(uniform_random_biv_poly(params_glwe, pol, precision),
	           "pol's malloc failed in uniform_random_biv_poly_dft.");

	biv_coefs_to_dft(module, params_glwe, result_dft, pol);

	status = 0;
cleanup:
	free(pol);
	return status;
}

//! COMMON PART (begin)

uint64_t poly_biv_bytes(const GLWEParams* params_glwe) { return poly_biv_coef_number(params_glwe) * sizeof(int64_t); }

uint64_t glwe_params_l(const GLWEParams* params_glwe) { return params_glwe->n_limbs / (params_glwe->k + 1); }

void biv_to_univ_rnx(const GLWEParams* params_glwe, double* res_univ, const PolyBiv* pol_biv)
{
	uint64_t nn      = params_glwe->nn;
	uint64_t kappa   = params_glwe->kappa;
	uint64_t l       = glwe_params_l(params_glwe);
	uint64_t l_max   = INT_ROUND_UP_DIV(53ul, kappa);
	uint64_t start_l = l > l_max ? l_max : l;

	// res_univ(X^p) = Sum_i{1,l}[poly(X^p, Y^i) * 2^(-kappa*i)]
	double pkappa = exp2(-(double)kappa);

	for (uint64_t p = 0; p < nn; p++)
	{
		res_univ[p] = 0;
		for (uint64_t i = start_l; i >= 1; --i)
		{
			res_univ[p] += (double)pol_biv[(i - 1) * nn + p];
			res_univ[p] *= pkappa;
		}
	}
}

int univ_rnx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivRnX* pol_univ)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t nn    = params_glwe->nn;
	uint64_t kappa = params_glwe->kappa;
	uint64_t l     = glwe_params_l(params_glwe);
	uint64_t l_max = INT_ROUND_UP_DIV(53ul, kappa);

	// Fills each pol_biv(X^p, Y^i) with the pol_univ's decomposition coefficients of  in [-2^(kappa* - 1) ; 2^(kappa -
	// 1) - 1]
	int64_t mask = (1LL << kappa) - 1;
	double acc   = 0;
	for (uint64_t i = 1; i <= l && i <= l_max; i++)
	{
		acc += ldexp(1.0, kappa - 1 - kappa * i);
	}

	for (uint64_t p = 0; p < nn; p++)
	{
		double tmp = pol_univ[p] + acc;

		if (tmp < 0) tmp -= floor(tmp);
		for (uint64_t i = 1; i <= l && i <= l_max; i++)
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
	uint64_t l  = glwe_params_l(params_glwe);
	pvda_vec_znx_dft(module, res_dft, l, a, l, nn);
	return 1;
}

int biv_dft_to_coefs(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBivDFT* a_dft)
{
	uint64_t nn = params_glwe->nn;
	uint64_t l  = glwe_params_l(params_glwe);
	return pvda_vec_znx_idft(module, res, l, a_dft, l);
}

int biv_to_univ_tnx(const GLWEParams* params_glwe, PolyUnivTnX* res_tnx, const PolyBiv* pol)
{
	uint64_t nn    = params_glwe->nn;
	uint64_t kappa = params_glwe->kappa;
	uint64_t l     = glwe_params_l(params_glwe);
	uint64_t l_max = INT_ROUND_UP_DIV(64ul, kappa);

	memset(res_tnx, 0, nn * sizeof(*res_tnx));

	for (uint64_t i = 0; i < l && i < l_max; ++i)
	{
		for (uint64_t p = 0; p < nn; p++)
		{
			int shft_amt     = 64 - (int)kappa - (int)(i * kappa);
			uint64_t add_amt = shft_amt > 0 ? ((uint64_t)pol[i * nn + p]) << shft_amt : (pol[i * nn + p] >> -shft_amt);
			res_tnx[p] += add_amt;
		}
	}
	return 0;
}

int univ_tnx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivTnX* pol_tnx)
{
	uint64_t nn    = params_glwe->nn;
	uint64_t kappa = params_glwe->kappa;
	uint64_t l     = glwe_params_l(params_glwe);
	uint64_t l_max = INT_ROUND_UP_DIV(64ul, kappa);

	uint64_t acc = 0;
	int64_t mask = (1LL << kappa) - 1;
	for (uint64_t i = 1; i <= l && i <= l_max; i++)
	{
		acc += 1LL << (kappa - 1 - kappa * i);
	}

	for (uint64_t p = 0; p < nn; p++)
	{
		uint64_t tmp = pol_tnx[p] + acc;

		for (uint64_t i = 1; i <= l && i <= l_max; i++)
		{
			int shft_amt = 64 - (int)(i * kappa);
			if (shft_amt > 0)
				res[(i - 1) * nn + p] = (int64_t)((tmp >> shft_amt) & mask) - (1LL << (kappa - 1));
			else
				res[(i - 1) * nn + p] = (int64_t)((tmp << -shft_amt) & mask) - (1LL << (kappa - 1));
		}
	}

	return 0;
}
