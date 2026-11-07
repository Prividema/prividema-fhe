#include "univariate_polynomial.h"

#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "backend.h"
#include "backend_arithmetic.h"

// =============================================
// |                                           |
// |                PolyUnivZnx                |
// |                                           |
// =============================================

PolyUnivZnx* new_univ_znx(const GLWEParams* params_glwe)
{
	// aligned_alloc guarantees the address to be a multiple of 64. This alignment matters for some backends.
	// ! WARNING : In C11 univ_znx_bytes(.) has to be a multiple of 64. Newer standard relaxed this.
	// ! WARNING : On Windows we use _aligned_malloc(size, 64) and free with _aligned free.
	return aligned_alloc(64, univ_znx_bytes(params_glwe));
}

void delete_univ_znx(PolyUnivZnx* pol) { free(pol); }

// =============================================
// |                                           |
// |                    DFT                    |
// |                                           |
// =============================================

BE_PolyUnivZnxDFT* new_univ_znx_dft(const PvdaBackend* module)
{
	return pvda_new_univ_znx_dft(module);
}

void delete_univ_znx_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* pol)
{
	pvda_delete_univ_znx_dft(module, pol);
}

int univ_znx_to_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* res_dft, const PolyUnivZnx* a)
{
	uint64_t nn  = pvda_module_extract_nn(module);
	FlatPolyVec a_vec = flat_poly_vec_view(nn, 1, nn, (int64_t*)a);
	return pvda_flat_univ_to_dft(module, res_dft, &a_vec);
}

int univ_dft_to_znx(const PvdaBackend* module, PolyUnivZnx* res, const BE_PolyUnivZnxDFT* a_dft)
{
	uint64_t nn    = pvda_module_extract_nn(module);
	FlatPolyVec res_vec = flat_poly_vec_view(nn, 1, nn, res);
	// TODO : "Unflat" the polynomial in res
	return pvda_flat_univ_idft(module, &res_vec, a_dft);
}

// =============================================
// |                                           |
// |                   Prep                    |
// |                                           |
// =============================================

BE_PolyUnivSvpPrep* new_univ_svp_prep(const PvdaBackend* module)
{
	return pvda_new_univ_svp_prep(module);
}

void delete_univ_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* pol)
{
	pvda_delete_univ_svp_prep(module, pol);
}

int univ_znx_to_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* res_prep, const PolyUnivZnx *a)
{
	uint64_t nn  = pvda_module_extract_nn(module);
	FlatPolyVec a_flat = flat_poly_vec_view(nn, 1, nn, (int64_t*)a);
	return pvda_flat_univ_to_svp_prep(module, res_prep, &a_flat);
}

// =============================================
// |                                           |
// |                   Tnx                     |
// |                                           |
// =============================================

PolyUnivTnX* new_univ_tnx(const GLWEParams* params_glwe) { return malloc(univ_tnx_bytes(params_glwe)); }

void delete_univ_tnx(PolyUnivTnX* pol) { free(pol); }

// =============================================
// |                                           |
// |                   Rnx                     |
// |                                           |
// =============================================

PolyUnivRnX* new_univ_rnx(const GLWEParams* params_glwe) { return malloc(univ_rnx_bytes(params_glwe)); }

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

uint64_t univ_znx_bytes(const GLWEParams* params_glwe) { return params_glwe->nn * sizeof(PolyUnivZnx);}

uint64_t univ_rnx_bytes(const GLWEParams* params_glwe) { return params_glwe->nn * sizeof(PolyUnivRnX); }

uint64_t univ_tnx_bytes(const GLWEParams* params_glwe) { return params_glwe->nn * sizeof(PolyUnivTnX); }

// =============================================
// |                                           |
// |              Print Functions              |
// |                                           |
// =============================================

void printf_poly_univ_ZnX(PolyUnivZnx* pol, int64_t nn)
{
	printf("%" PRId64, pol[0]);
	for (int64_t p = 1; p < nn; p++)
	{
		if (pol[p] < 0) printf(" - %" PRIu64 " X^%" PRId64, -(uint64_t)pol[p], p);
		else            printf(" + %" PRId64 " X^%" PRId64, pol[p], p);
	}
	putchar('\n');
}

void printf_poly_univ_RnX(double* pol, int64_t nn)
{
	printf("%g", pol[0]);
	for (int64_t p = 1; p < nn; p++)
	{
		if (pol[p] < 0) printf(" - %g X^%" PRId64, -pol[p], p);
		else            printf(" + %g X^%" PRId64, pol[p], p);
	}
	putchar('\n');
}

void printf_vec_poly_univ_TnX(int64_t* pol, int64_t nn)
{
	printf("%" PRIu64, pol[0]);
	for (uint64_t p = 1; p < nn; p++)
		printf(" + %" PRIu64 " X^%" PRIu64, pol[p], p);
	putchar('\n');
}
