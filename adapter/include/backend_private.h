#ifndef PVDA_BACKEND_PRIVATE_H
#define PVDA_BACKEND_PRIVATE_H

#include <stdint.h>

#include "backend.h"
#include "backend_arithmetic.h"
#include "backend_nt.h"
#include "backend_rng.h"
#include "pvda_ffts.h"
#include "pvda_types.h"

/**
 * @file backend_private.h
 *
 * @brief Internals of the adapter: the virtual tables that dispatch each adapter function to the implementation
 * of the selected backend, and the backend object itself.
 *
 * A backend holds two virtual tables:
 * - vt: polynomial arithmetic (spqlios, cpu number theory, ...);
 * - rng_vt: random number generation (cpu, openssl), chosen independently of vt.
 */

typedef struct module_info_t SPQLIOS_MODULE;

typedef typeof(pvda_new_univ_znx_dft) PVDA_NEW_UNIV_ZNX_DFT_F;
typedef typeof(pvda_delete_univ_znx_dft) PVDA_DELETE_UNIV_ZNX_DFT_F;
typedef typeof(pvda_univ_znx_to_dft) PVDA_UNIV_ZNX_TO_DFT_F;
typedef typeof(pvda_univ_znx_idft) PVDA_UNIV_ZNX_IDFT_F;

typedef typeof(pvda_new_univ_svp_prep) PVDA_NEW_UNIV_SVP_PREP_F;
typedef typeof(pvda_delete_univ_svp_prep) PVDA_DELETE_UNIV_SVP_PREP_F;
typedef typeof(pvda_univ_znx_to_svp_prep) PVDA_UNIV_ZNX_TO_SVP_PREP_F;

typedef typeof(pvda_new_biv_dft) PVDA_NEW_BIV_DFT_F;
typedef typeof(pvda_delete_biv_dft) PVDA_DELETE_BIV_DFT_F;
typedef typeof(pvda_biv_to_dft) PVDA_BIV_TO_DFT_F;
typedef typeof(pvda_biv_idft) PVDA_BIV_IDFT_F;
typedef typeof(pvda_biv_dft_at_const) PVDA_BIV_DFT_AT_F;
typedef typeof(pvda_univ_dft_as_biv_dft_const) PVDA_UNIV_DFT_AS_BIV_DFT_F;

typedef typeof(pvda_svp_apply_dft) PVDA_SVP_APPLY_DFT_F;
typedef typeof(pvda_svp_apply_dft_to_dft) PVDA_SVP_APPLY_DFT_TO_DFT_F;

typedef typeof(pvda_new_biv_vmp_prep) PVDA_NEW_BIV_VMP_PREP_F;
typedef typeof(pvda_delete_biv_vmp_prep) PVDA_DELETE_BIV_VMP_PREP_F;
typedef typeof(pvda_biv_to_vmp_prep) PVDA_BIV_TO_VMP_PREP_F;

typedef typeof(pvda_new_vmp_pmat) PVDA_NEW_VMP_PMAT_F;
typedef typeof(pvda_biv_to_vmp_pmat) PVDA_BIV_TO_VMP_PMAT_F;
typedef typeof(pvda_delete_vmp_pmat) PVDA_DELETE_VMP_PMAT_F;

typedef typeof(pvda_vmp_apply_dft) PVDA_VMP_APPLY_DFT_F;
typedef typeof(pvda_vmp_apply_dft_to_dft) PVDA_VMP_APPLY_DFT_TO_DFT_F;
typedef typeof(pvda_vmp_apply_prepared_to_dft) PVDA_VMP_APPLY_PREPARED_TO_DFT_F;

typedef typeof(pvda_biv_normalize_base2k) PVDA_BIV_NORMALIZE_BASE2K_F;
typedef typeof(pvda_univ_znx_small_product) PVDA_UNIV_ZNX_SMALL_PRODUCT_F;
typedef typeof(pvda_biv_negate) PVDA_BIV_NEGATE_F;
typedef typeof(pvda_biv_add) PVDA_BIV_ADD_F;
typedef typeof(pvda_biv_sub) PVDA_BIV_SUB_F;
typedef typeof(pvda_univ_znx_automorphism) PVDA_UNIV_ZNX_AUTOMORPHISM_F;
typedef typeof(pvda_biv_automorphism) PVDA_BIV_AUTOMORPHISM_F;
typedef typeof(pvda_biv_rotate) PVDA_BIV_ROTATE_F;
typedef typeof(pvda_module_extract_nn) PVDA_MODULE_EXTRACT_NN_F;
typedef typeof(pvda_ring_2nth_root) PVDA_RING_2NTH_ROOT_F;
typedef typeof(pvda_biv_dft_add) PVDA_BIV_DFT_ADD_F;

typedef struct pvda_virtual_table
{
	PVDA_NEW_UNIV_ZNX_DFT_F* pvda_new_univ_znx_dft;
	PVDA_DELETE_UNIV_ZNX_DFT_F* pvda_delete_univ_znx_dft;
	PVDA_UNIV_ZNX_TO_DFT_F* pvda_univ_znx_to_dft;
	PVDA_UNIV_ZNX_IDFT_F* pvda_univ_znx_idft;
	PVDA_UNIV_DFT_AS_BIV_DFT_F* pvda_univ_dft_as_biv_dft;

	PVDA_NEW_UNIV_SVP_PREP_F* pvda_new_univ_svp_prep;
	PVDA_DELETE_UNIV_SVP_PREP_F* pvda_delete_univ_svp_prep;
	PVDA_UNIV_ZNX_TO_SVP_PREP_F* pvda_univ_znx_to_svp_prep;

	PVDA_NEW_BIV_DFT_F* pvda_new_biv_dft;
	PVDA_DELETE_BIV_DFT_F* pvda_delete_biv_dft;
	PVDA_BIV_TO_DFT_F* pvda_biv_to_dft;
	PVDA_BIV_IDFT_F* pvda_biv_idft;
	PVDA_BIV_DFT_AT_F* pvda_biv_dft_at;

	PVDA_SVP_APPLY_DFT_F* pvda_svp_apply_dft;
	PVDA_SVP_APPLY_DFT_TO_DFT_F* pvda_svp_apply_dft_to_dft;

	PVDA_NEW_VMP_PMAT_F* pvda_new_vmp_pmat;
	PVDA_BIV_TO_VMP_PMAT_F* pvda_biv_to_vmp_pmat;
	PVDA_DELETE_VMP_PMAT_F* pvda_delete_vmp_pmat;

	PVDA_NEW_BIV_VMP_PREP_F* pvda_new_biv_vmp_prep;
	PVDA_DELETE_BIV_VMP_PREP_F* pvda_delete_biv_vmp_prep;
	PVDA_BIV_TO_VMP_PREP_F* pvda_biv_to_vmp_prep;

	PVDA_VMP_APPLY_DFT_F* pvda_vmp_apply_dft;
	PVDA_VMP_APPLY_DFT_TO_DFT_F* pvda_vmp_apply_dft_to_dft;
	PVDA_VMP_APPLY_PREPARED_TO_DFT_F* pvda_vmp_apply_prepared_to_dft;

	PVDA_BIV_NORMALIZE_BASE2K_F* pvda_biv_normalize_base2k;
	PVDA_UNIV_ZNX_SMALL_PRODUCT_F* pvda_univ_znx_small_product;
	PVDA_BIV_NEGATE_F* pvda_biv_negate;
	PVDA_BIV_ADD_F* pvda_biv_add;
	PVDA_BIV_SUB_F* pvda_biv_sub;
	PVDA_UNIV_ZNX_AUTOMORPHISM_F* pvda_univ_znx_automorphism;
	PVDA_BIV_AUTOMORPHISM_F* pvda_biv_automorphism;
	PVDA_BIV_ROTATE_F* pvda_biv_rotate;
	PVDA_MODULE_EXTRACT_NN_F* pvda_module_extract_nn;
	PVDA_RING_2NTH_ROOT_F* pvda_ring_2nth_root;
	PVDA_BIV_DFT_ADD_F *pvda_biv_dft_add;
} PvdaVirtualTable;

typedef typeof(pvda_rng_bytes) PVDA_RNG_BYTES_F;
typedef typeof(pvda_rng_uniform_pow2_vec) PVDA_RNG_UNIFORM_POW2_VEC_F;
typedef typeof(pvda_rng_uniform_range_vec) PVDA_RNG_UNIFORM_RANGE_VEC_F;
typedef typeof(pvda_rng_binary_vec) PVDA_RNG_BINARY_VEC_F;
typedef typeof(pvda_rng_centered_binomial_vec) PVDA_RNG_CENTERED_BINOMIAL_VEC_F;
typedef typeof(pvda_rng_normal_vec) PVDA_RNG_NORMAL_VEC_F;

/**
 * @brief RNG virtual table.
 *
 * An RNG implementation must provide pvda_rng_bytes. The distribution entries are first filled with the generic
 * implementations of rng_generic.c (built on pvda_rng_bytes); an implementation can then override any of them.
 */
typedef struct pvda_rng_virtual_table
{
	PVDA_RNG_BYTES_F* pvda_rng_bytes;  ///< mandatory
	PVDA_RNG_UNIFORM_POW2_VEC_F* pvda_rng_uniform_pow2_vec;
	PVDA_RNG_UNIFORM_RANGE_VEC_F* pvda_rng_uniform_range_vec;
	PVDA_RNG_BINARY_VEC_F* pvda_rng_binary_vec;
	PVDA_RNG_CENTERED_BINOMIAL_VEC_F* pvda_rng_centered_binomial_vec;
	PVDA_RNG_NORMAL_VEC_F* pvda_rng_normal_vec;

	/** @brief Frees the RNG state (NULL if the implementation has no state). */
	void (*delete_state)(void* state);
} PvdaRngVirtualTable;

struct pvda_module_t
{
	SPQLIOS_MODULE* spqlios_module;
	struct pvda_fft_data_t* pvda_fft_data;

	PvdaVirtualTable vt;  ///< Arithmetic

	PvdaRngVirtualTable rng_vt;  ///< Randomness
	PvdaRngType rng_type;        ///< RNG implementation behind rng_vt
	int rng_seeded;              ///< 1 if deterministic (seeded), 0 if drawing from the system source
	void* rng_state;             ///< Implementation-specific state (e.g. keystream state when seeded), owned
};

/** @brief Fills the distribution entries of vt with the generic implementations built on vt->pvda_rng_bytes. */
void pvda_rng_fill_generic(PvdaRngVirtualTable* vt);

/** @brief Frees the RNG state of a backend and clears its RNG virtual table. */
void pvda_delete_backend_rng(PvdaBackend* module);

#endif
