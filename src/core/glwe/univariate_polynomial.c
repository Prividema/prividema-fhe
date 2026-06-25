#include "univariate_polynomial.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "backend.h"
#include "bivariate_polynomial.h"

int univ_coefs_to_dft(const PVDA_MODULE* module, PolyUnivDFT* res_dft, const PolyUniv* a)
{
	uint64_t nn   = pvda_module_extract_nn(module);
	PolyBiv a_biv = new_biv_view(nn, 1, nn, a);
	pvda_vec_znx_dft(module, res_dft, 1, &a_biv);
	return 0;
}

int univ_dft_to_coefs(const PVDA_MODULE* module, PolyUniv* res, const PolyUnivDFT* a_dft)
{
	uint64_t nn     = pvda_module_extract_nn(module);
	PolyBiv res_biv = new_biv_view(nn, 1, nn, res);
	return pvda_vec_znx_idft(module, &res_biv, a_dft, 1);
}

uint64_t poly_univ_bytes(const GLWEParams* params_glwe)
{
	uint64_t nn = params_glwe->nn;
	return nn * sizeof(PolyUniv);
}

uint64_t poly_univ_rnx_bytes(const GLWEParams* params_glwe) { return params_glwe->nn * sizeof(PolyUnivRnX); }

uint64_t poly_univ_tnx_bytes(const GLWEParams* params_glwe) { return params_glwe->nn * sizeof(PolyUnivTnX); }

PolyUniv* new_univ(const GLWEParams* params_glwe) { return aligned_alloc(64, poly_univ_bytes(params_glwe)); }

void delete_univ(PolyUniv* pol) { free(pol); }

PolyUnivRnX* new_univ_rnx(const GLWEParams* params_glwe) { return malloc(poly_univ_rnx_bytes(params_glwe)); }

void delete_univ_rnx(PolyUnivRnX* pol) { free(pol); }

PolyUnivDFT* new_univ_dft(const PVDA_MODULE* module) { return pvda_new_vec_znx_dft(module, 1); }

void delete_univ_dft(const PVDA_MODULE* module, PolyUnivDFT* pol) { pvda_delete_vec_znx_dft(module, pol); /*TODO fix*/ }

PolyUnivTnX* new_univ_tnx(const GLWEParams* params_glwe) { return malloc(poly_univ_tnx_bytes(params_glwe)); }

void delete_univ_tnx(PolyUnivTnX* pol) { free(pol); }

int univ_rnx_to_tnx(const GLWEParams* params_glwe, PolyUnivTnX* res, PolyUnivRnX* a)
{
	for (size_t i = 0; i < params_glwe->nn; ++i)
	{
		res[i] = (uint64_t)(ldexp(a[i] - floor(a[i]), 64));
	}
	return 0;
}

int univ_tnx_to_rnx(const GLWEParams* params_glwe, PolyUnivRnX* res, PolyUnivTnX* a)
{
	for (size_t i = 0; i < params_glwe->nn; ++i)
	{
		res[i] = ldexp((double)((int64_t)a[i]), -64);
	}
	return 0;
}
