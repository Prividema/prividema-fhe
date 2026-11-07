#ifndef PVDA_BACKEND_ARITH
#define PVDA_BACKEND_ARITH

#include "backend.h"
#include "backend_types.h"

/**
 * @file backend_arithmetic.h
 * 
 * @brief Arithmetic backend functions for prividema.
 * 
 * This file contains all the arithmetic functions from the backend we need on the library side.
 */

// =============================================
// |                                           |
// |              BE_PolyUnivZnxDFT            |
// |                                           |
// =============================================

BE_PolyUnivZnxDFT* pvda_new_univ_znx_dft(const PvdaBackend* module);
void pvda_delete_univ_znx_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* pol);
int pvda_flat_univ_to_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* res, const FlatPolyVec* a);
int pvda_flat_univ_idft(const PvdaBackend* module, FlatPolyVec* res, const BE_PolyUnivZnxDFT* pol);

// =============================================
// |                                           |
// |             BE_PolyUnivSvpPrep            |
// |                                           |
// =============================================

BE_PolyUnivSvpPrep* pvda_new_univ_svp_prep(const PvdaBackend* module);
void pvda_delete_univ_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* pol);
int pvda_flat_univ_to_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* res, const FlatPolyVec *a);

// =============================================
// |                                           |
// |                BE_PolyBivDFT              |
// |                                           |
// =============================================

BE_PolyBivDFT* pvda_new_biv_dft(const PvdaBackend* module, uint64_t size);
void pvda_delete_biv_dft(const PvdaBackend* module, BE_PolyBivDFT* res);

int pvda_flat_biv_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size, const FlatPolyVec* a);
int pvda_flat_biv_idft(const PvdaBackend* module, FlatPolyVec* res, const BE_PolyBivDFT* pol, uint64_t pol_size);

// =============================================
// |                                           |
// |                    SVP                    |
// |                                           |
// =============================================

int pvda_svp_apply_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size,
                       const BE_PolyUnivSvpPrep* prepared_pol, const FlatPolyVec* a);

int pvda_svp_apply_dft_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size,
                              const BE_PolyUnivSvpPrep* ppol, const BE_PolyBivDFT* a, uint64_t a_size);

//ZnxBig
int64_t* pvda_new_vec_znx_big(const PvdaBackend* module, int64_t size);
void pvda_delete_vec_znx_big(const PvdaBackend* module, int64_t* res);

// =============================================
// |                                           |
// |             BE_PolyBivVmpPrep             |
// |                                           |
// =============================================

BE_PolyBivVmpPrep* pvda_new_biv_vmp_prep(const PvdaBackend* module, uint64_t biv_l);
void pvda_delete_biv_vmp_prep(const PvdaBackend* module, BE_PolyBivVmpPrep* biv_prep);
int pvda_flat_biv_to_vmp_prep(const PvdaBackend* module, BE_PolyBivVmpPrep* res, uint64_t res_size, const FlatPolyVec* a);

// =============================================
// |                                           |
// |              BE_PolyMatVmpPrep            |
// |                                           |
// =============================================

BE_PolyMatVmpPrep* pvda_new_vmp_pmat(const PvdaBackend* module, uint64_t nrows, uint64_t ncols);
int pvda_vmp_prepare_contiguous(const PvdaBackend* module, BE_PolyMatVmpPrep* pmat, const int64_t* mat, uint64_t nrows, uint64_t ncols);
void pvda_delete_vmp_pmat(const PvdaBackend* module, BE_PolyMatVmpPrep* pmat);

// =============================================
// |                                           |
// |                    VMP                    |
// |                                           |
// =============================================

int pvda_vmp_apply_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size, const FlatPolyVec* a, const BE_PolyMatVmpPrep* pmat, uint64_t nrows, uint64_t ncols);
int pvda_vmp_apply_dft_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, const uint64_t res_size,
                              const BE_PolyBivDFT* a_dft, uint64_t a_size,
                              const BE_PolyMatVmpPrep* pmat, const uint64_t nrows, const uint64_t ncols);
int pvda_vmp_apply_prepared_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, const uint64_t res_size,
                                   const BE_PolyBivVmpPrep* a_prep, uint64_t a_size,
                                   const BE_PolyMatVmpPrep* pmat, const uint64_t nrows, const uint64_t ncols);

// =============================================
// |                                           |
// |                 Arithmetic                |
// |                                           |
// =============================================

int pvda_vec_znx_normalize_base2k(const PvdaBackend* module, uint64_t log2_base2k, FlatPolyVec* res, const FlatPolyVec* a);
int pvda_znx_small_product(const PvdaBackend* module, int64_t* res, const int64_t* a, const int64_t* b);
int pvda_vec_znx_negate(const PvdaBackend* module, FlatPolyVec* res, const FlatPolyVec* a);
int pvda_vec_znx_add(const PvdaBackend* module, FlatPolyVec* res, const FlatPolyVec* a, const FlatPolyVec* b);
int pvda_vec_znx_sub(const PvdaBackend* module, FlatPolyVec* res, const FlatPolyVec* a, const FlatPolyVec* b);
int pvda_znx_automorphism(const PvdaBackend* module, const int64_t p, int64_t* res, const int64_t* a);
int pvda_vec_znx_automorphism(const PvdaBackend* module, const int64_t p, FlatPolyVec* res, const FlatPolyVec* a);
int pvda_vec_znx_rotate(const PvdaBackend* module, const int64_t p, FlatPolyVec* res, const FlatPolyVec* a);
int pvda_vec_znx_dft_add(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_n_limbs, BE_PolyBivDFT* a,
                         uint64_t a_n_limbs, BE_PolyBivDFT* b, uint64_t b_n_limbs);
#endif
