#ifndef SPQLIOS_ALIAS_H
#define SPQLIOS_ALIAS_H

#include <stdint.h>

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

//typedef enum module_type_t MODULE_TYPE;

/// @brief Opaque structure that describe the modules (\f$\mathbb{Z}_n[X]\f$, \f$\mathbb{T}_n[X]\f$) and the hardware.
typedef struct module_info_t MODULE;

/// @brief Opaque type that represents a prepared matrix.
typedef struct vmp_pmat_t VMP_PMAT;

/// @brief Opaque type that represents a vector of \f$\mathbb{Z}_n[X]\f$ in the DFT domain.
typedef struct vec_znx_dft_t VEC_ZNX_DFT;

/// @brief Opaque type that represents a vector of \f$\mathbb{Z}_n[X]\f$ in large coeffs space.
typedef struct vec_znx_bigcoeff_t VEC_ZNX_BIG;

/// @brief Opaque type that represents a prepared scalar vector product.
typedef struct svp_ppol_t SVP_PPOL;

/// @brief Opaque type that represents a prepared left convolution vector product.
typedef struct cnv_pvec_l_t CNV_PVEC_L;

/// @brief Opaque type that represents a prepared right convolution vector product.
typedef struct cnv_pvec_r_t CNV_PVEC_R;

// =============================================
// |                                           |
// |             spqlios Functions             |
// |                                           |
// |       The following function are all      |
// |      named pvda_f where f is the name of  |
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
MODULE* pvda_new_module_info(uint64_t nn);

void pvda_delete_module_info(MODULE* module);

VecUnivDFT* pvda_new_vec_znx_dft(const MODULE* module, uint64_t size);

void pvda_vec_znx_dft(const MODULE* module, double* res, uint64_t res_size, const int64_t* a, uint64_t a_size,
                      uint64_t a_sl);

void pvda_delete_vec_znx_dft(double* res);

int64_t* pvda_new_vec_znx_big(const MODULE* module, int64_t size);

int pvda_vec_znx_idft(const MODULE* module, int64_t* res, uint64_t res_size, const double* a_dft, uint64_t a_size);

void pvda_delete_vec_znx_big(int64_t* res);

double* pvda_new_svp_ppol(const MODULE* module);

void pvda_svp_prepare(const MODULE* module, PolyUnivDFT* prepared_pol, const int64_t* pol);

void pvda_svp_apply_dft(const MODULE* module, const double* res, uint64_t res_size, const PolyUnivDFT* prepared_pol,
                        const int64_t* a, uint64_t a_size, uint64_t a_sl);

void pvda_delete_svp_ppol(double* res);

double* pvda_new_vmp_pmat(const MODULE* module, uint64_t nrows, uint64_t ncols);

int pvda_vmp_prepare_contiguous(const MODULE* module, double* pmat, const int64_t* mat, uint64_t nrows, uint64_t ncols);

int pvda_vmp_apply_dft(const MODULE* module, double* res, uint64_t res_size, const int64_t* a, uint64_t a_size,
                       uint64_t a_sl, const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols);

int pvda_vmp_apply_dft_to_dft(const MODULE* module, VecBivDFT* res, const uint64_t res_size, const VecBivDFT* a_dft,
                              uint64_t a_size, const MatBivDFT* pmat, const uint64_t nrows, const uint64_t ncols);

void pvda_delete_vmp_pmat(double* pmat);

int pvda_vec_znx_normalize_base2k(const MODULE* module, uint64_t log2_base2k, int64_t* res, uint64_t res_size,
                                  uint64_t res_sl, const int64_t* a, uint64_t a_size, uint64_t a_sl);

int pvda_znx_small_product(const MODULE* module, int64_t* res, const int64_t* a, const int64_t* b);

int pvda_vec_znx_negate(const MODULE* module, int64_t* res, uint64_t res_size, uint64_t res_sl, const int64_t* a,
                        uint64_t a_size, uint64_t a_sl);

int pvda_vec_znx_add(const MODULE* module, int64_t* res, uint64_t res_size, uint64_t res_sl, const int64_t* a,
                     uint64_t a_size, uint64_t a_sl, const int64_t* b, uint64_t b_size, uint64_t b_sl);

int pvda_vec_znx_sub(const MODULE* module, int64_t* res, uint64_t res_size, uint64_t res_sl, const int64_t* a,
                     uint64_t a_size, uint64_t a_sl, const int64_t* b, uint64_t b_size, uint64_t b_sl);

int pvda_vec_znx_automorphism(const MODULE* module, const int64_t p, int64_t* res, uint64_t res_size, uint64_t res_sl,
                              const int64_t* a, uint64_t a_size, uint64_t a_sl);

uint64_t pvda_module_extract_nn(const MODULE* module);

/**@}*/
#endif
