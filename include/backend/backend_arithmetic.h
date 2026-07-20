#ifndef PVDA_BACKEND_ARITH
#define PVDA_BACKEND_ARITH

#include "backend.h"
#include "maths_structures.h"

VecUnivDFT* pvda_new_vec_znx_dft(const PvdaBackend* backend, uint64_t size);

int pvda_vec_znx_dft(const PvdaBackend* backend, PolyBivDFT* res, uint64_t res_size, const PolyBiv* a);

void pvda_delete_vec_znx_dft(const PvdaBackend* backend, PolyBivDFT* res);

int64_t* pvda_new_vec_znx_big(const PvdaBackend* backend, int64_t size);

int pvda_vec_znx_idft(const PvdaBackend* backend, PolyBiv* res, const PolyBivDFT* a_dft, uint64_t a_size);

void pvda_delete_vec_znx_big(const PvdaBackend* backend, int64_t* res);

PolyUnivDFT* pvda_new_svp_ppol(const PvdaBackend* backend);

int pvda_svp_prepare(const PvdaBackend* backend, PolyUnivDFT* prepared_pol, const int64_t* pol);

int pvda_svp_apply_dft(const PvdaBackend* backend, const PolyBivDFT* res, uint64_t res_size,
                       const PolyUnivDFT* prepared_pol, const PolyBiv* a);

int pvda_svp_apply_dft_to_dft(const PvdaBackend* backend, const PolyBivDFT* res, uint64_t res_size,
                              const PolyUnivDFT* ppol, const PolyBivDFT* a, uint64_t a_size);

void pvda_delete_svp_ppol(const PvdaBackend* backend, PolyUnivDFT* res);

MatBivDFT* pvda_new_vmp_pmat(const PvdaBackend* backend, uint64_t nrows, uint64_t ncols);

int pvda_vmp_prepare_contiguous(const PvdaBackend* backend, MatBivDFT* pmat, const int64_t* mat, uint64_t nrows,
                                uint64_t ncols);

int pvda_vmp_apply_dft(const PvdaBackend* backend, PolyBivDFT* res, uint64_t res_size, const PolyBiv* a,
                       const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols);

int pvda_vmp_apply_dft_to_dft(const PvdaBackend* backend, VecBivDFT* res, const uint64_t res_size,
                              const VecBivDFT* a_dft, uint64_t a_size, const MatBivDFT* pmat, const uint64_t nrows,
                              const uint64_t ncols);

int pvda_vmp_apply_prepared_to_dft(const PvdaBackend* backend, VecBivDFT* res, const uint64_t res_size,
                                   const VecBivDFT* a_dft, uint64_t a_size, const MatBivDFT* pmat, const uint64_t nrows,
                                   const uint64_t ncols);

void pvda_delete_vmp_pmat(const PvdaBackend* backend, double* pmat);

int pvda_vec_znx_normalize_base2k(const PvdaBackend* backend, uint64_t log2_base2k, PolyBiv* res, const PolyBiv* a);

int pvda_znx_small_product(const PvdaBackend* backend, PolyUniv* res, const PolyUniv* a, const PolyUniv* b);

int pvda_vec_znx_negate(const PvdaBackend* backend, PolyBiv* res, const PolyBiv* a);

int pvda_vec_znx_add(const PvdaBackend* backend, PolyBiv* res, const PolyBiv* a, const PolyBiv* b);

int pvda_vec_znx_sub(const PvdaBackend* backend, PolyBiv* res, const PolyBiv* a, const PolyBiv* b);

int pvda_znx_automorphism(const PvdaBackend* backend, const int64_t p, PolyUniv* res, const PolyUniv* a);
int pvda_vec_znx_automorphism(const PvdaBackend* backend, const int64_t p, PolyBiv* res, const PolyBiv* a);

int pvda_vec_znx_rotate(const PvdaBackend* backend, const int64_t p, PolyBiv* res, const PolyBiv* a);

int pvda_vmp_prepare_vec(const PvdaBackend* backend, PolyBivPrep* pvec, uint64_t nrows, const PolyBiv* a);
#endif
