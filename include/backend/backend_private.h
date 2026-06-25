#ifndef PVDA_BACKEND_PRIVATE_H
#define PVDA_BACKEND_PRIVATE_H

#include <stdint.h>

#include "backend.h"

typedef struct module_info_t SPQLIOS_MODULE;

typedef typeof(pvda_vmp_apply_dft) PVDA_VMP_APPLY_DFT_F;
typedef typeof(pvda_new_vec_znx_dft) PVDA_NEW_VEC_ZNX_DFT_F;
typedef typeof(pvda_vec_znx_dft) PVDA_VEC_ZNX_DFT_F;
typedef typeof(pvda_delete_vec_znx_dft) PVDA_DELETE_VEC_ZNX_DFT_F;
typedef typeof(pvda_new_vec_znx_big) PVDA_NEW_VEC_ZNX_BIG_F;
typedef typeof(pvda_vec_znx_idft) PVDA_VEC_ZNX_IDFT_F;
typedef typeof(pvda_delete_vec_znx_big) PVDA_DELETE_VEC_ZNX_BIG_F;
typedef typeof(pvda_new_svp_ppol) PVDA_NEW_SVP_PPOL_F;
typedef typeof(pvda_svp_prepare) PVDA_SVP_PREPARE_F;
typedef typeof(pvda_svp_apply_dft) PVDA_SVP_APPLY_DFT_F;
typedef typeof(pvda_svp_apply_dft_to_dft) PVDA_SVP_APPLY_DFT_TO_DFT_F;
typedef typeof(pvda_delete_svp_ppol) PVDA_DELETE_SVP_PPOL_F;
typedef typeof(pvda_new_vmp_pmat) PVDA_NEW_VMP_PMAT_F;
typedef typeof(pvda_vmp_prepare_contiguous) PVDA_VMP_PREPARE_CONTIGUOUS_F;
typedef typeof(pvda_vmp_apply_dft) PVDA_VMP_APPLY_DFT_F;
typedef typeof(pvda_vmp_apply_dft_to_dft) PVDA_VMP_APPLY_DFT_TO_DFT_F;
typedef typeof(pvda_vmp_apply_prepared_to_dft) PVDA_VMP_APPLY_PREPARED_TO_DFT_F;
typedef typeof(pvda_delete_vmp_pmat) PVDA_DELETE_VMP_PMAT_F;
typedef typeof(pvda_vec_znx_normalize_base2k) PVDA_VEC_ZNX_NORMALIZE_BASE2K_F;
typedef typeof(pvda_znx_small_product) PVDA_ZNX_SMALL_PRODUCT_F;
typedef typeof(pvda_vec_znx_negate) PVDA_VEC_ZNX_NEGATE_F;
typedef typeof(pvda_vec_znx_add) PVDA_VEC_ZNX_ADD_F;
typedef typeof(pvda_vec_znx_sub) PVDA_VEC_ZNX_SUB_F;
typedef typeof(pvda_znx_automorphism) PVDA_ZNX_AUTOMORPHISM_F;
typedef typeof(pvda_vec_znx_automorphism) PVDA_VEC_ZNX_AUTOMORPHISM_F;
typedef typeof(pvda_vec_znx_rotate) PVDA_VEC_ZNX_ROTATE_F;
typedef typeof(pvda_module_extract_nn) PVDA_MODULE_EXTRACT_NN_F;
typedef typeof(pvda_vmp_prepare_vec) PVDA_VMP_PREPARE_VEC_F;

typedef struct pvda_virtual_table
{
	PVDA_NEW_VEC_ZNX_DFT_F* pvda_new_vec_znx_dft;
	PVDA_VEC_ZNX_DFT_F* pvda_vec_znx_dft;
	PVDA_DELETE_VEC_ZNX_DFT_F* pvda_delete_vec_znx_dft;
	PVDA_NEW_VEC_ZNX_BIG_F* pvda_new_vec_znx_big;
	PVDA_VEC_ZNX_IDFT_F* pvda_vec_znx_idft;
	PVDA_DELETE_VEC_ZNX_BIG_F* pvda_delete_vec_znx_big;
	PVDA_NEW_SVP_PPOL_F* pvda_new_svp_ppol;
	PVDA_SVP_PREPARE_F* pvda_svp_prepare;
	PVDA_SVP_APPLY_DFT_F* pvda_svp_apply_dft;
	PVDA_SVP_APPLY_DFT_TO_DFT_F* pvda_svp_apply_dft_to_dft;
	PVDA_DELETE_SVP_PPOL_F* pvda_delete_svp_ppol;
	PVDA_NEW_VMP_PMAT_F* pvda_new_vmp_pmat;
	PVDA_VMP_PREPARE_CONTIGUOUS_F* pvda_vmp_prepare_contiguous;
	PVDA_VMP_APPLY_DFT_F* pvda_vmp_apply_dft;
	PVDA_VMP_APPLY_DFT_TO_DFT_F* pvda_vmp_apply_dft_to_dft;
	PVDA_VMP_APPLY_PREPARED_TO_DFT_F* pvda_vmp_apply_prepared_to_dft;
	PVDA_DELETE_VMP_PMAT_F* pvda_delete_vmp_pmat;
	PVDA_VEC_ZNX_NORMALIZE_BASE2K_F* pvda_vec_znx_normalize_base2k;
	PVDA_ZNX_SMALL_PRODUCT_F* pvda_znx_small_product;
	PVDA_VEC_ZNX_NEGATE_F* pvda_vec_znx_negate;
	PVDA_VEC_ZNX_ADD_F* pvda_vec_znx_add;
	PVDA_VEC_ZNX_SUB_F* pvda_vec_znx_sub;
	PVDA_ZNX_AUTOMORPHISM_F* pvda_znx_automorphism;
	PVDA_VEC_ZNX_AUTOMORPHISM_F* pvda_vec_znx_automorphism;
	PVDA_VEC_ZNX_ROTATE_F* pvda_vec_znx_rotate;
	PVDA_MODULE_EXTRACT_NN_F* pvda_module_extract_nn;
	PVDA_VMP_PREPARE_VEC_F* pvda_vmp_prepare_vec;
} PvdaVirutalTable;

struct pvda_module_t
{
	SPQLIOS_MODULE* spqlios_module;

	PvdaVirutalTable vt;
};

typedef struct pvda_module_t PvdaBackend;

#endif
