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
 * @brief Internals of the adapter: the virtual table that dispatches each adapter function to the implementation
 * of the selected backend, and the backend object itself.
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
typedef typeof(pvda_rand_uniform_pow2) PVDA_RAND_UNIFORM_POW2_F;
typedef typeof(pvda_rand_uniform_pow2_vec) PVDA_RAND_UNIFORM_POW2_VEC_F;
typedef typeof(pvda_rand_uniform_binary_vec) PVDA_RAND_UNIFORM_BINARY_VEC_F;
typedef typeof(pvda_rand_uniform) PVDA_RAND_UNIFORM_F;
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
	PVDA_RAND_UNIFORM_POW2_F* pvda_rand_uniform_pow2;
	PVDA_RAND_UNIFORM_POW2_VEC_F* pvda_rand_uniform_pow2_vec;
	PVDA_RAND_UNIFORM_BINARY_VEC_F* pvda_rand_uniform_binary_vec;
	PVDA_RAND_UNIFORM_F* pvda_rand_uniform;
	PVDA_RING_2NTH_ROOT_F* pvda_ring_2nth_root;
	PVDA_BIV_DFT_ADD_F *pvda_biv_dft_add;
} PvdaVirtualTable;

struct pvda_module_t
{
	SPQLIOS_MODULE* spqlios_module;
	struct pvda_fft_data_t* pvda_fft_data;

	PvdaVirtualTable vt;
};

typedef struct pvda_module_t PvdaBackend;

#endif
