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
	// C11 requires the size to be a multiple of the alignment: round it up (it matters when nn < 8).
	// ! WARNING : On Windows we use _aligned_malloc(size, 64) and free with _aligned free.
	uint64_t bytes = univ_znx_bytes(params_glwe);
	return aligned_alloc(64, (bytes + 63) / 64 * 64);
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
	return pvda_univ_znx_to_dft(module, res_dft, a);
}

int univ_dft_to_znx(const PvdaBackend* module, PolyUnivZnx* res, const BE_PolyUnivZnxDFT* a_dft)
{
	return pvda_univ_znx_idft(module, res, a_dft);
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
	return pvda_univ_znx_to_svp_prep(module, res_prep, a);
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
// |               PolyUnivZnxMP               |
// |                                           |
// =============================================

PolyUnivZnxMP* new_univ_znx_mp(uint64_t nn, uint64_t nb_words)
{
	if (nn == 0 || nb_words == 0) return NULL;

	PolyUnivZnxMP* pol = malloc(sizeof(PolyUnivZnxMP));
	if (!pol) return NULL;

	pol->nn       = nn;
	pol->nb_words = nb_words;
	pol->ptr      = calloc(nn * nb_words, sizeof(uint64_t));
	if (!pol->ptr)
	{
		free(pol);
		return NULL;
	}
	return pol;
}

void delete_univ_znx_mp(PolyUnivZnxMP* pol)
{
	if (!pol) return;
	free(pol->ptr);
	free(pol);
}

uint64_t univ_znx_mp_bytes(const PolyUnivZnxMP* pol) { return pol->nn * pol->nb_words * sizeof(uint64_t); }

int univ_znx_to_znx_mp(PolyUnivZnxMP* res, const PolyUnivZnx* a)
{
	for (uint64_t p = 0; p < res->nn; ++p)
	{
		uint64_t* coef = res->ptr + p * res->nb_words;
		uint64_t ext   = a[p] < 0 ? UINT64_MAX : 0;
		coef[0]        = (uint64_t)a[p];
		for (uint64_t w = 1; w < res->nb_words; ++w) coef[w] = ext;
	}
	return 0;
}

int univ_znx_mp_to_znx(PolyUnivZnx* res, const PolyUnivZnxMP* a)
{
	int status = 0;
	for (uint64_t p = 0; p < a->nn; ++p)
	{
		const uint64_t* coef = a->ptr + p * a->nb_words;
		res[p]               = (int64_t)coef[0];

		// The value fits in an int64_t iff every upper word is the sign extension of word 0
		uint64_t ext = res[p] < 0 ? UINT64_MAX : 0;
		for (uint64_t w = 1; w < a->nb_words; ++w)
			if (coef[w] != ext) status = -1;
	}
	return status;
}

// =============================================
// |                                           |
// |            Tnx-Rnx Conversions            |
// |                                           |
// =============================================

int univ_rnx_to_tnx(const GLWEParams* params_glwe, PolyUnivTnX* res, PolyUnivRnX* a)
{
	for (size_t i = 0; i < params_glwe->nn; ++i)
	{
		// a - floor(a) is in [0, 1], but rounds to exactly 1 for tiny negative a: 2^64 does not fit in a uint64_t
		// (the conversion would be undefined behaviour) and stands for 0 on the torus.
		double scaled = ldexp(a[i] - floor(a[i]), 64);
		res[i]        = scaled < 0x1p64 ? (uint64_t)scaled : 0;
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

void print_poly_univ_ZnX(PolyUnivZnx* pol, int64_t nn)
{
	printf("%" PRId64, pol[0]);
	for (int64_t p = 1; p < nn; p++)
	{
		if (pol[p] < 0) printf(" - %" PRIu64 " X^%" PRId64, -(uint64_t)pol[p], p);
		else            printf(" + %" PRId64 " X^%" PRId64, pol[p], p);
	}
	putchar('\n');
}

void print_poly_univ_RnX(PolyUnivRnX* pol, int64_t nn)
{
	printf("%g", pol[0]);
	for (int64_t p = 1; p < nn; p++)
	{
		if (pol[p] < 0) printf(" - %g X^%" PRId64, -pol[p], p);
		else            printf(" + %g X^%" PRId64, pol[p], p);
	}
	putchar('\n');
}

void print_poly_univ_TnX(PolyUnivTnX* pol, int64_t nn)
{
	printf("%" PRIu64, pol[0]);
	for (int64_t p = 1; p < nn; p++)
		printf(" + %" PRIu64 " X^%" PRId64, pol[p], p);
	putchar('\n');
}

void print_poly_univ_ZnX_mp(const PolyUnivZnxMP* pol)
{
	for (uint64_t p = 0; p < pol->nn; p++)
	{
		printf("X^%" PRIu64 ": 0x", p);
		for (uint64_t w = pol->nb_words; w-- > 0;) printf("%016" PRIx64, pol->ptr[p * pol->nb_words + w]);
		putchar('\n');
	}
}
