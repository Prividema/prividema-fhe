#ifndef PVDA_BACKEND_ARITH
#define PVDA_BACKEND_ARITH

#include "backend.h"
#include "backend_types.h"
#include "maths_structures.h"

/**
 * @file backend_arithmetic.h
 * @brief Arithmetic backend functions for prividema
 *
 * Refer to spqlios' documentation for the functions here
 */

BE_PolyVecDFT* pvda_new_vec_znx_dft(const PvdaBackend* backend, uint64_t size);

int pvda_vec_znx_dft(const PvdaBackend* backend, BE_PolyVecDFT* res, uint64_t res_size, const ZnxVec* a);

void pvda_delete_vec_znx_dft(const PvdaBackend* backend, BE_PolyVecDFT* res);

int64_t* pvda_new_vec_znx_big(const PvdaBackend* backend, int64_t size);

int pvda_vec_znx_idft(const PvdaBackend* backend, ZnxVec* res, const BE_PolyVecDFT* a_dft, uint64_t a_size);

void pvda_delete_vec_znx_big(const PvdaBackend* backend, int64_t* res);

BE_PolySvpPrep* pvda_new_svp_ppol(const PvdaBackend* backend);

int pvda_svp_prepare(const PvdaBackend* backend, BE_PolySvpPrep* prepared_pol, const int64_t* pol);

int pvda_svp_apply_dft(const PvdaBackend* backend, const BE_PolyVecDFT* res, uint64_t res_size,
                       const BE_PolySvpPrep* prepared_pol, const ZnxVec* a);

int pvda_svp_apply_dft_to_dft(const PvdaBackend* backend, const BE_PolyVecDFT* res, uint64_t res_size,
                              const BE_PolySvpPrep* ppol, const BE_PolyVecDFT* a, uint64_t a_size);

void pvda_delete_svp_ppol(const PvdaBackend* backend, BE_PolySvpPrep* ppol);

BE_PolyVecVmpPrep* pvda_new_vmp_pvec(const PvdaBackend* backend, uint64_t nrows);

void pvda_delete_vmp_pvec(const PvdaBackend* backend, BE_PolyVecVmpPrep* pvec);

BE_PolyMatVmpPrep* pvda_new_vmp_pmat(const PvdaBackend* backend, uint64_t nrows, uint64_t ncols);

int pvda_vmp_prepare_contiguous(const PvdaBackend* backend, BE_PolyMatVmpPrep* pmat, const int64_t* mat, uint64_t nrows,
                                uint64_t ncols);

int pvda_vmp_apply_dft(const PvdaBackend* backend, BE_PolyVecDFT* res, uint64_t res_size, const ZnxVec* a,
                       const BE_PolyMatVmpPrep* pmat, uint64_t nrows, uint64_t ncols);

int pvda_vmp_apply_dft_to_dft(const PvdaBackend* backend, BE_PolyVecDFT* res, const uint64_t res_size,
                              const BE_PolyVecDFT* a_dft, uint64_t a_size, const BE_PolyMatVmpPrep* pmat, const uint64_t nrows,
                              const uint64_t ncols);

int pvda_vmp_apply_prepared_to_dft(const PvdaBackend* backend, BE_PolyVecDFT* res, const uint64_t res_size,
                                   const BE_PolyVecVmpPrep* a_prep, uint64_t a_size, const BE_PolyMatVmpPrep* pmat, const uint64_t nrows,
                                   const uint64_t ncols);

void pvda_delete_vmp_pmat(const PvdaBackend* backend, BE_PolyMatVmpPrep* pmat);

int pvda_vec_znx_normalize_base2k(const PvdaBackend* backend, uint64_t log2_base2k, ZnxVec* res, const ZnxVec* a);

int pvda_znx_small_product(const PvdaBackend* backend, PolyUniv* res, const PolyUniv* a, const PolyUniv* b);

int pvda_vec_znx_negate(const PvdaBackend* backend, ZnxVec* res, const ZnxVec* a);

int pvda_vec_znx_add(const PvdaBackend* backend, ZnxVec* res, const ZnxVec* a, const ZnxVec* b);

int pvda_vec_znx_sub(const PvdaBackend* backend, ZnxVec* res, const ZnxVec* a, const ZnxVec* b);

int pvda_znx_automorphism(const PvdaBackend* backend, const int64_t p, PolyUniv* res, const PolyUniv* a);
int pvda_vec_znx_automorphism(const PvdaBackend* backend, const int64_t p, ZnxVec* res, const ZnxVec* a);

int pvda_vec_znx_rotate(const PvdaBackend* backend, const int64_t p, ZnxVec* res, const ZnxVec* a);

int pvda_vmp_prepare_vec(const PvdaBackend* backend, BE_PolyVecVmpPrep* pvec, uint64_t nrows, const ZnxVec* a);

int pvda_vec_znx_dft_add(const PvdaBackend* backend, BE_PolyVecDFT* res, uint64_t res_n_limbs,
                         BE_PolyVecDFT* a, uint64_t a_n_limbs,
                         BE_PolyVecDFT* b, uint64_t b_n_limbs);
#endif
