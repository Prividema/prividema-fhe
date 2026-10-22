#include "univariate_polynomial.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "backend.h"
#include "backend_arithmetic.h"

// =============================================
// |                                           |
// |                 PolyUniv                  |
// |                                           |
// =============================================

PolyUniv* new_univ(const GLWEParams* params_glwe)
{
	return aligned_alloc(64, univ_bytes(params_glwe));
}

void delete_univ(PolyUniv* pol) { free(pol); }

// =============================================
// |                                           |
// |                    DFT                    |
// |                                           |
// =============================================

BE_PolyVecDFT* new_univ_dft(const PvdaBackend* module)
{
	return pvda_new_vec_znx_dft(module, 1);
}

void delete_univ_dft(const PvdaBackend* module, BE_PolyVecDFT* pol)
{
	pvda_delete_vec_znx_dft(module, pol);
}

int univ_coefs_to_dft(const PvdaBackend* module, BE_PolyVecDFT* res_dft, const PolyUniv* a)
{
	uint64_t nn  = pvda_module_extract_nn(module);
	ZnxVec a_vec = znx_vec_view(nn, 1, nn, (int64_t*)a);
	pvda_vec_znx_dft(module, res_dft, 1, &a_vec);
	return 0;
}

int univ_dft_to_coefs(const PvdaBackend* module, PolyUniv* res, const BE_PolyVecDFT* a_dft)
{
	uint64_t nn    = pvda_module_extract_nn(module);
	ZnxVec res_vec = znx_vec_view(nn, 1, nn, res);
	return pvda_vec_znx_idft(module, &res_vec, a_dft, 1);
}

// =============================================
// |                                           |
// |                   Prep                    |
// |                                           |
// =============================================

BE_PolySvpPrep* new_univ_prep(const PvdaBackend* module)
{
	return pvda_new_svp_ppol(module);
}

void delete_univ_prep(const PvdaBackend* module, BE_PolySvpPrep* pol)
{
	pvda_delete_svp_ppol(module, pol);
}

int univ_coefs_to_prep(const PvdaBackend* module, BE_PolySvpPrep* res_prep, const PolyUniv *a)
{
	return pvda_svp_prepare(module, res_prep, a);
}

// =============================================
// |                                           |
// |                   Tnx                     |
// |                                           |
// =============================================

PolyUnivTnX* new_univ_tnx(const GLWEParams* params_glwe)
{
	return malloc(univ_tnx_bytes(params_glwe));
}

void delete_univ_tnx(PolyUnivTnX* pol) { free(pol); }

// =============================================
// |                                           |
// |                   Rnx                     |
// |                                           |
// =============================================

PolyUnivRnX* new_univ_rnx(const GLWEParams* params_glwe)
{
	return malloc(univ_rnx_bytes(params_glwe));
}

void delete_univ_rnx(PolyUnivRnX* pol) { free(pol); }

// =============================================
// |                                           |
// |            Tnx-Rnx Conversions            |
// |                                           |
// =============================================

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

// =============================================
// |                                           |
// |                   Utils                   |
// |                                           |
// =============================================

uint64_t univ_bytes(const GLWEParams* params_glwe)
{
	uint64_t nn = params_glwe->nn;
	return nn * sizeof(PolyUniv);
}

uint64_t univ_rnx_bytes(const GLWEParams* params_glwe) { return params_glwe->nn * sizeof(PolyUnivRnX); }

uint64_t univ_tnx_bytes(const GLWEParams* params_glwe) { return params_glwe->nn * sizeof(PolyUnivTnX); }

// =============================================
// |                                           |
// |              Print Functions              |
// |                                           |
// =============================================

void printf_poly_univ_ZnX(PolyUniv* pol, int64_t nn)
{
	printf((pol[0] < 0) ? "%ld " : " %ld ", pol[0]);
	for (uint64_t p = 1; p < nn; p++)
		printf((pol[p] < 0) ? "- %ld X^%ld " : "+ %ld X^%ld ", (pol[p] < 0) ? -pol[p] : pol[p], p);
}

void printf_poly_univ_RnX(double* pol, int64_t nn)
{
	for (uint64_t p = 0; p < nn; p++) printf(" %lf X^%ld ", pol[p], p);
}

void printf_vec_poly_univ(int64_t* pols, int64_t pols_size, int64_t nn)
{
	printf("\n");
	for (uint64_t j = 0; j < pols_size; j++)
	{
		printf("\n%ld-th component : ", j);
		printf_poly_univ_ZnX(pols + j * nn, nn);
		printf("\n");
	}
}
