#include "univariate_polynomial.h"

#include <math.h>
#include <stdint.h>

#include "rng.h"

int univ_coefs_to_dft(const MODULE* module, PolyUnivDFT* res_dft, const PolyUniv* a)
{
	pvda_vec_znx_dft(module, res_dft, 1, a, 1, pvda_module_extract_nn(module));
	return 1;
}

int univ_dft_to_coefs(const MODULE* module, PolyUniv* res, const PolyUnivDFT* a_dft)
{
	return pvda_vec_znx_idft(module, res, 1, a_dft, pvda_module_extract_nn(module));
}

uint64_t poly_univ_bytes(const GLWEParams* params_glwe)
{
	uint64_t nn = params_glwe->nn;
	return nn * sizeof(PolyUniv);
}

uint64_t poly_univ_rnx_bytes(const GLWEParams* params_glwe) { return params_glwe->nn * sizeof(PolyUnivRnX); }

uint64_t poly_univ_tnx_bytes(const GLWEParams* params_glwe) { return params_glwe->nn * sizeof(PolyUnivTnX); }

PolyUniv* new_univ(const GLWEParams* params_glwe) { return malloc(poly_univ_bytes(params_glwe)); }

void delete_univ(PolyUniv* pol) { free(pol); }

PolyUnivRnX* new_univ_rnx(const GLWEParams* params_glwe) { return malloc(poly_univ_rnx_bytes(params_glwe)); }

void delete_univ_rnx(PolyUnivRnX* pol) { free(pol); }

PolyUnivDFT* new_univ_dft(const MODULE* module) { return pvda_new_vec_znx_dft(module, 1); }

void delete_univ_dft(PolyUnivDFT* pol) { pvda_delete_vec_znx_dft(pol); }

PolyUnivTnX* new_univ_tnx(const GLWEParams* params_glwe) { return malloc(poly_univ_tnx_bytes(params_glwe)); }

void delete_univ_tnx(PolyUnivTnX* pol) { free(pol); }

int univ_rnx_to_tnx(const GLWEParams* params_glwe, PolyUnivTnX* res, PolyUnivRnX* a)
{
	for (size_t i = 0; i < params_glwe->nn; ++i)
	{
		res[i] = (uint64_t)(ldexp(a[i] - floor(a[i]), 64));
	}
	return 1;
}

int univ_tnx_to_rnx(const GLWEParams* params_glwe, PolyUnivRnX* res, PolyUnivTnX* a)
{
	for (size_t i = 0; i < params_glwe->nn; ++i)
	{
		res[i] = ldexp((double)((int64_t)a[i]), -64);
	}
	return 1;
}
