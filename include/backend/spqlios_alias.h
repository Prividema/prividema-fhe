#ifndef SPQLIOS_ALIAS_H
#define SPQLIOS_ALIAS_H

#include <stdint.h>

#include "backend.h"
#include "backend_private.h"
#include "maths_structures.h"

/**
 *
 * @file spqlios_alias.h
 *
 *
 * Provides a very thin wrapper over spqlios functions
 *
 *
 * @note For documentation on the functions in this file, refer to the splqios library
 *
 */

// =============================================
// |                                           |
// |      Aliases for spqlios structures       |
// |                                           |
// =============================================

/// @brief Opaque type that represents a prepared matrix.
typedef struct vmp_pmat_t VMP_PMAT;

/// @brief Opaque type that represents a vector of \f$\mathbb{Z}_n[X]\f$ in the DFT domain.
typedef struct vec_znx_dft_t VEC_ZNX_DFT;

/// @brief Opaque type that represents a vector of \f$\mathbb{Z}_n[X]\f$ in large coeffs space.
typedef struct vec_znx_bigcoeff_t VEC_ZNX_BIG;

/// @brief Opaque type that represents a prepared scalar vector product.
typedef struct svp_ppol_t SVP_PPOL;

// =============================================
// |                                           |
// |             spqlios Functions             |
// |                                           |
// |       The following function are all      |
// |      named spqlios_f where f is the name of  |
// |           a function as it might exist    |
// |      in spqlios. Sometimes it does not    |
// |      exist and we have a temporary        |
// |                implementation.            |
// |                                           |
// |  We additionnally allocate the memory of  |
// |                the structures.            |
// |                                           |
// =============================================
/**@{*/
/** @name Spqlios wrapper functions (see spqlios docs) */
SPQLIOS_MODULE* spqlios_new_module_info(uint64_t nn);

void spqlios_delete_module_info(PvdaModule* module);

VecUnivDFT* spqlios_new_vec_znx_dft(const PvdaModule* module, uint64_t size);

int spqlios_vec_znx_dft(const PvdaModule* module, double* res, uint64_t res_size, const PolyBiv* a);

void spqlios_delete_vec_znx_dft(const PvdaModule* module, double* res);

int64_t* spqlios_new_vec_znx_big(const PvdaModule* module, int64_t size);

int spqlios_vec_znx_idft(const PvdaModule* module, PolyBiv* res, const double* a_dft, uint64_t a_size);

void spqlios_delete_vec_znx_big(const PvdaModule* module, int64_t* res);

double* spqlios_new_svp_ppol(const PvdaModule* module);

int spqlios_svp_prepare(const PvdaModule* module, PolyUnivDFT* prepared_pol, const int64_t* pol);

int spqlios_svp_apply_dft(const PvdaModule* module, const double* res, uint64_t res_size,
                          const PolyUnivDFT* prepared_pol, const PolyBiv* a);

int spqlios_svp_apply_dft_to_dft(const PvdaModule* module, const double* res, uint64_t res_size,
                                 const PolyUnivDFT* ppol, const PolyBivDFT* a, uint64_t a_size);

void spqlios_delete_svp_ppol(const PvdaModule* module, double* res);

double* spqlios_new_vmp_pmat(const PvdaModule* module, uint64_t nrows, uint64_t ncols);

int spqlios_vmp_prepare_contiguous(const PvdaModule* module, double* pmat, const int64_t* mat, uint64_t nrows,
                                   uint64_t ncols);

int spqlios_vmp_apply_dft(const PvdaModule* module, double* res, uint64_t res_size, const PolyBiv* a,
                          const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols);

int spqlios_vmp_apply_dft_to_dft(const PvdaModule* module, VecBivDFT* res, const uint64_t res_size,
                                 const VecBivDFT* a_dft, uint64_t a_size, const MatBivDFT* pmat, const uint64_t nrows,
                                 const uint64_t ncols);

int spqlios_vmp_apply_prepared_to_dft(const PvdaModule* module, VecBivDFT* res, const uint64_t res_size,
                                      const VecBivDFT* a_dft, uint64_t a_size, const MatBivDFT* pmat,
                                      const uint64_t nrows, const uint64_t ncols);

void spqlios_delete_vmp_pmat(const PvdaModule* module, double* pmat);

int spqlios_vec_znx_normalize_base2k(const PvdaModule* module, uint64_t log2_base2k, PolyBiv* res, const PolyBiv* a);

int spqlios_znx_small_product(const PvdaModule* module, PolyUniv* res, const PolyUniv* a, const PolyUniv* b);

int spqlios_vec_znx_negate(const PvdaModule* module, PolyBiv* res, const PolyBiv* a);

int spqlios_vec_znx_add(const PvdaModule* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b);

int spqlios_vec_znx_sub(const PvdaModule* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b);

int spqlios_znx_automorphism(const PvdaModule* module, const int64_t p, PolyUniv* res, const PolyUniv* a);
int spqlios_vec_znx_automorphism(const PvdaModule* module, const int64_t p, PolyBiv* res, const PolyBiv* a);

int spqlios_vec_znx_rotate(const PvdaModule* module, const int64_t p, PolyBiv* res, const PolyBiv* a);

uint64_t spqlios_module_extract_nn(const PvdaModule* module);

int spqlios_vmp_prepare_vec(const PvdaModule* module, double* pvec, uint64_t nrows, const PolyBiv* a);

/**@}*/
#endif
