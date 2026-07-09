#ifndef PVDA_BACKEND_H
#define PVDA_BACKEND_H

#include <stdint.h>
#include <sys/types.h>

#include "maths_structures.h"

typedef struct pvda_module_t PvdaBackend;

typedef struct pvda_module_config_t
{
	uint64_t nn;

} PvdaBackendConfig;

//Legacy backend port functions
PvdaBackend* pvda_new_spqlios_backend(int nn);

//Functions

/* SPQLIOS */

void pvda_delete_backend(PvdaBackend* backend);

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

uint64_t pvda_module_extract_nn(const PvdaBackend* backend);

int pvda_vmp_prepare_vec(const PvdaBackend* backend, PolyBivPrep* pvec, uint64_t nrows, const PolyBiv* a);

//RNG

/**
 * @brief Generates a uniformly sampled random number in [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param result  The resulting uniformly sampled integer
 * @param nb_bits The number of bits of the result
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int pvda_rand_uniform_pow2(const PvdaBackend* backend, int64_t* result, uint64_t nb_bits);

/**
 * @brief Generates a uniformly sampled vector of random numbers in [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param result  The resulting uniformly sampled integer
 * @param n       Number of elements in the vector
 * @param nb_bits The number of bits of the result
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int pvda_rand_uniform_pow2_vec(const PvdaBackend* backend, int64_t* result, uint64_t n, uint64_t nb_bits);

/**
 * @brief Generates a uniformly sampled vector of random numbers in {0,1}
 *
 * @param result  The resulting uniformly sampled integer
 * @param n       Number of elements in the vector
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int pvda_rand_uniform_binary_vec(const PvdaBackend* backend, uint64_t* result, uint64_t n);

/**
 * @brief Generates a uniformly sampled random number in [limit_down, limit_up] via
 * power-of-2 sampling and resampling if out-of-bounds
 *
 * @param module      The backend object
 * @param result      The resulting uniformly sampled integer
 * @param limit_down  The lower bound of the uniform sample
 * @param limit_up    The upper bound of the uniform sample
 *
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int pvda_rand_uniform(const PvdaBackend* module, int64_t* result, int64_t limit_down, int64_t limit_up);

#endif
