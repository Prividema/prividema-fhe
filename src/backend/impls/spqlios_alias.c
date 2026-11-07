#include "impls/spqlios_alias.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "backend.h"
#include "utils.h"
#include "vec_znx_arithmetic.h"

SPQLIOS_MODULE* spqlios_new_module_info(uint64_t nn) { return new_module_info(nn, FFT64); }

void spqlios_delete_module_info(PvdaBackend* module)
{
	if (!module || !module->spqlios_module) return;
	delete_module_info(module->spqlios_module);
}

// =============================================
// |                                           |
// |              BE_PolyUnivZnxDFT            |
// |                                           |
// =============================================

BE_PolyUnivZnxDFT* spqlios_new_univ_znx_dft(const PvdaBackend* module)
{
	return (BE_PolyUnivZnxDFT*)new_vec_znx_dft(module->spqlios_module, 1); 
}

void spqlios_delete_univ_znx_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* pol)
{
	if (!pol) return;
	delete_vec_znx_dft((VEC_ZNX_DFT*)pol);
}

int spqlios_flat_univ_to_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* res, const FlatPolyVec* a)
{
	vec_znx_dft(module->spqlios_module, (VEC_ZNX_DFT*)res, 1, a->ptr, a->size, a->stride);
	return 0;
}

int spqlios_flat_univ_idft(const PvdaBackend* module, FlatPolyVec* res, const BE_PolyUnivZnxDFT* pol)
{
	int status = -1;

	size_t tmp_bytes = vec_znx_idft_tmp_bytes(module->spqlios_module);
	uint8_t* tmp     = malloc(tmp_bytes);
	if (tmp_bytes) CHECK_ALLOC(tmp, "ntt120_vec_znx_idft_tmp_bytes_avx failed in vec_znx_idft_p");

	assert(res->stride == res->nn);
	vec_znx_idft(module->spqlios_module, (VEC_ZNX_BIG*)res->ptr, 1, (VEC_ZNX_DFT*)pol, 1, tmp);

	status = 0;

cleanup:
	free(tmp);
	return status;
}

// =============================================
// |                                           |
// |             BE_PolyUnivSvpPrep            |
// |                                           |
// =============================================

BE_PolyUnivSvpPrep* spqlios_new_univ_svp_prep(const PvdaBackend* module)
{
	return (BE_PolyUnivSvpPrep*)new_svp_ppol(module->spqlios_module);
}

void spqlios_delete_univ_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* ppol)
{
	if (!ppol) return;
	delete_svp_ppol(((SVP_PPOL*)ppol));
}

int spqlios_flat_univ_to_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* ppol, const FlatPolyVec* pol)
{
	svp_prepare(module->spqlios_module, (SVP_PPOL*)ppol, pol->ptr);
	return 0;
}

// =============================================
// |                                           |
// |                BE_PolyBivDFT              |
// |                                           |
// =============================================

BE_PolyBivDFT* spqlios_new_biv_dft(const PvdaBackend* module, uint64_t size)
{
	return (BE_PolyBivDFT*)new_vec_znx_dft(module->spqlios_module, size); 
}

void spqlios_delete_biv_dft(const PvdaBackend* module, BE_PolyBivDFT* pol)
{
	if (!pol) return;
	delete_vec_znx_dft((VEC_ZNX_DFT*)pol);
}

int spqlios_flat_biv_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size, const FlatPolyVec* a)
{
	vec_znx_dft(module->spqlios_module, (VEC_ZNX_DFT*)res, res_size, a->ptr, a->size, a->stride);
	return 0;
}

int spqlios_flat_biv_idft(const PvdaBackend* module, FlatPolyVec* res, const BE_PolyBivDFT* a_dft, uint64_t a_size)
{
	int status = -1;

	size_t tmp_bytes = vec_znx_idft_tmp_bytes(module->spqlios_module);
	uint8_t* tmp     = malloc(tmp_bytes);
	if (tmp_bytes) CHECK_ALLOC(tmp, "ntt120_vec_znx_idft_tmp_bytes_avx failed in vec_znx_idft_p");

	assert(res->stride == res->nn);
	vec_znx_idft(module->spqlios_module, (VEC_ZNX_BIG*)res->ptr, res->size, (VEC_ZNX_DFT*)a_dft, a_size, tmp);

	status = 0;

cleanup:
	free(tmp);
	return status;
}

// =============================================
// |                                           |
// |                    SVP                    |
// |                                           |
// =============================================

int spqlios_svp_apply_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size,
                          const BE_PolyUnivSvpPrep* prepared_pol, const FlatPolyVec* a)
{
	svp_apply_dft(module->spqlios_module, (VEC_ZNX_DFT*)res, res_size, (SVP_PPOL*)prepared_pol, a->ptr, a->size,
	              a->stride);
	return 0;
}

int spqlios_svp_apply_dft_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size,
                              const BE_PolyUnivSvpPrep* ppol, const BE_PolyBivDFT* a, uint64_t a_size)
{
	svp_apply_dft_to_dft(module->spqlios_module, (VEC_ZNX_DFT*)res, res_size, (SVP_PPOL*)ppol, (VEC_ZNX_DFT*)a, a_size);
	return 0;
}

// ZnxBig

int64_t* spqlios_new_vec_znx_big(const PvdaBackend* module, int64_t size)
{
	return (int64_t*)new_vec_znx_big(module->spqlios_module, size);
}

void spqlios_delete_vec_znx_big(const PvdaBackend* module, int64_t* res)
{
	if (!res) return;
	delete_vec_znx_big((VEC_ZNX_BIG*)res);
}

// =============================================
// |                                           |
// |             BE_PolyBivVmpPrep             |
// |                                           |
// =============================================

BE_PolyBivVmpPrep* spqlios_new_biv_vmp_prep(const PvdaBackend* module, uint64_t nrows)
{
	uint64_t bytes = module_get_n(module->spqlios_module) * nrows * sizeof(double);
	void* res = spqlios_alloc(bytes);
	if (res) memset(res, 0, bytes);
	return (BE_PolyBivVmpPrep*)res;
}

void spqlios_delete_biv_vmp_prep(const PvdaBackend* module, BE_PolyBivVmpPrep* pvec)
{
	if (!pvec) return;
	spqlios_free(pvec);
}

int spqlios_flat_biv_to_vmp_prep(const PvdaBackend* module, BE_PolyBivVmpPrep* pvec, uint64_t nrows, const FlatPolyVec* a)
{
	int status = -1;

	void* tmp_space = aligned_alloc(64, vmp_prepare_vector_tmp_bytes(module->spqlios_module, nrows, a->size));
	CHECK_ALLOC(tmp_space, "tmp space alloc failed in spqlios_convolution_prepare");

	vmp_prepare_vector(module->spqlios_module, (VMP_PVEC*)pvec, nrows, a->ptr, a->size, a->stride, tmp_space);
	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

// =============================================
// |                                           |
// |              BE_PolyMatVmpPrep            |
// |                                           |
// =============================================

BE_PolyMatVmpPrep* spqlios_new_vmp_pmat(const PvdaBackend* module, uint64_t nrows, uint64_t ncols)
{
	return (BE_PolyMatVmpPrep*)new_vmp_pmat(module->spqlios_module, nrows, ncols);
}

int spqlios_vmp_prepare_contiguous(const PvdaBackend* module, BE_PolyMatVmpPrep* pmat, const int64_t* mat, uint64_t nrows,
                                   uint64_t ncols)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vmp_prepare_contiguous_tmp_bytes(module->spqlios_module, nrows, ncols));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vmp_prepare_contiguous_p");

	vmp_prepare_contiguous(module->spqlios_module, (VMP_PMAT*)pmat, mat, nrows, ncols, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

void spqlios_delete_vmp_pmat(const PvdaBackend* module, BE_PolyMatVmpPrep* pmat)
{
	if (!pmat) return;
	delete_vmp_pmat(((VMP_PMAT*)pmat));
}

// =============================================
// |                                           |
// |                    VMP                    |
// |                                           |
// =============================================
int spqlios_vmp_apply_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size, const FlatPolyVec* a, const BE_PolyMatVmpPrep* pmat, uint64_t nrows, uint64_t ncols)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vmp_apply_dft_tmp_bytes(module->spqlios_module, res_size, a->size, nrows, ncols));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vmp_apply_dft_p");

	vmp_apply_dft(module->spqlios_module, (VEC_ZNX_DFT*)res, res_size, a->ptr, a->size, a->stride, (VMP_PMAT*)pmat,
	              nrows, ncols, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

int spqlios_vmp_apply_dft_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, const uint64_t res_size,
                              	 const BE_PolyBivDFT* a_dft, uint64_t a_size, const BE_PolyMatVmpPrep* pmat,
                              	 const uint64_t nrows, const uint64_t ncols)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vmp_apply_dft_to_dft_tmp_bytes(module->spqlios_module, res_size, a_size, nrows, ncols));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vmp_apply_dft_to_dft_p");

	vmp_apply_dft_to_dft(module->spqlios_module, (VEC_ZNX_DFT*)res, res_size, (VEC_ZNX_DFT*)a_dft, a_size,
	                     (VMP_PMAT*)pmat, nrows, ncols, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

int spqlios_vmp_apply_prepared_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, const uint64_t res_size,
                                   	  const BE_PolyBivVmpPrep* a_prep, uint64_t a_size,
                                   	  const BE_PolyMatVmpPrep* pmat, const uint64_t nrows, const uint64_t ncols)
{
	int status = -1;
	uint8_t* tmp_space =
	    malloc(vmp_apply_pvec_to_dft_tmp_bytes(module->spqlios_module, res_size, a_size, nrows, ncols));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vmp_apply_dft_to_dft_p");

	vmp_apply_pvec_to_dft(module->spqlios_module, (VEC_ZNX_DFT*)res, res_size, (VMP_PVEC*)a_prep, a_size,
	                      (VMP_PMAT*)pmat, nrows, ncols, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

// =============================================
// |                                           |
// |                 Arithmetic                |
// |                                           |
// =============================================

// spqlios declares vec_dft_add but does not implement it, so it is implemented here.
// In the FFT64 layout, each DFT polynomial is nn doubles and the DFT is linear,
// so the addition is coefficient-wise.
// Same size convention as vec_znx_add: add up to the smallest size, copy up to the largest, then zero-extend.
int spqlios_vec_znx_dft_add(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_n_limbs, BE_PolyBivDFT* a,
                            uint64_t a_n_limbs, BE_PolyBivDFT* b, uint64_t b_n_limbs)
{
	uint64_t nn = module_get_n(module->spqlios_module);

	double* r       = (double*)res;
	const double* x = (const double*)a;
	const double* y = (const double*)b;

	// Make x the shortest input
	if (a_n_limbs > b_n_limbs)
	{
		const double* tmp = x;
		x                 = y;
		y                 = tmp;
		uint64_t tmp_n    = a_n_limbs;
		a_n_limbs         = b_n_limbs;
		b_n_limbs         = tmp_n;
	}

	uint64_t sum_n  = res_n_limbs < a_n_limbs ? res_n_limbs : a_n_limbs;
	uint64_t copy_n = res_n_limbs < b_n_limbs ? res_n_limbs : b_n_limbs;

	for (uint64_t i = 0; i < sum_n * nn; ++i) r[i] = x[i] + y[i];
	if (copy_n > sum_n) memmove(r + sum_n * nn, y + sum_n * nn, (copy_n - sum_n) * nn * sizeof(double));
	if (res_n_limbs > copy_n) memset(r + copy_n * nn, 0, (res_n_limbs - copy_n) * nn * sizeof(double));

	return 0;
}

int spqlios_vec_znx_normalize_base2k(const PvdaBackend* module, uint64_t log2_base2k, FlatPolyVec* res, const FlatPolyVec* a)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vec_znx_normalize_base2k_tmp_bytes(module->spqlios_module));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vec_znx_normalize_base2k_p");

	vec_znx_normalize_base2k(module->spqlios_module, log2_base2k, res->ptr, res->size, res->stride, a->ptr, a->size,
	                         a->stride, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

int spqlios_znx_small_product(const PvdaBackend* module, int64_t* res, const int64_t* a, const int64_t* b)
{
	int status         = -1;
	size_t tmp_size    = znx_small_single_product_tmp_bytes(module->spqlios_module);
	uint8_t* tmp_space = malloc(tmp_size);
	if (tmp_size) CHECK_ALLOC(tmp_space, "failed temp space alloc for polynomial multiplication");
	znx_small_single_product(module->spqlios_module, res, a, b, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

int spqlios_vec_znx_negate(const PvdaBackend* module, FlatPolyVec* res, const FlatPolyVec* a)
{
	vec_znx_negate(module->spqlios_module, res->ptr, res->size, res->stride, a->ptr, a->size, a->stride);
	return 0;
}

uint64_t spqlios_module_extract_nn(const PvdaBackend* module) { return module_get_n(module->spqlios_module); }

int spqlios_vec_znx_add(const PvdaBackend* module, FlatPolyVec* res, const FlatPolyVec* a, const FlatPolyVec* b)
{
	vec_znx_add(module->spqlios_module, res->ptr, res->size, res->stride, a->ptr, a->size, a->stride, b->ptr, b->size,
	            b->stride);
	return 0;
}
int spqlios_vec_znx_sub(const PvdaBackend* module, FlatPolyVec* res, const FlatPolyVec* a, const FlatPolyVec* b)
{
	vec_znx_sub(module->spqlios_module, res->ptr, res->size, res->stride, a->ptr, a->size, a->stride, b->ptr, b->size,
	            b->stride);
	return 0;
}

int spqlios_znx_automorphism(const PvdaBackend* module, const int64_t p, int64_t* res, const int64_t* a)
{
	vec_znx_automorphism(module->spqlios_module, p, res, 1, 0, a, 1, 0);
	return 0;
}
int spqlios_vec_znx_automorphism(const PvdaBackend* module, const int64_t p, FlatPolyVec* res, const FlatPolyVec* a)
{
	vec_znx_automorphism(module->spqlios_module, p, res->ptr, res->size, res->stride, a->ptr, a->size, a->stride);
	return 0;
}

int spqlios_vec_znx_rotate(const PvdaBackend* module, const int64_t p, FlatPolyVec* res, const FlatPolyVec* a)
{
	vec_znx_rotate(module->spqlios_module, p, res->ptr, res->size, res->stride, a->ptr, a->size, a->stride);
	return 0;
}

void pvda_fill_spqlios(struct pvda_virtual_table* vt)
{
	vt->pvda_new_univ_znx_dft    = spqlios_new_univ_znx_dft;
	vt->pvda_delete_univ_znx_dft = spqlios_delete_univ_znx_dft;
	vt->pvda_flat_univ_to_dft    = spqlios_flat_univ_to_dft;
	vt->pvda_flat_univ_idft      = spqlios_flat_univ_idft;

	vt->pvda_new_univ_svp_prep     = spqlios_new_univ_svp_prep;
	vt->pvda_delete_univ_svp_prep  = spqlios_delete_univ_svp_prep;
	vt->pvda_flat_univ_to_svp_prep = spqlios_flat_univ_to_svp_prep;

	vt->pvda_new_biv_dft     = spqlios_new_biv_dft;
	vt->pvda_delete_biv_dft  = spqlios_delete_biv_dft;
	vt->pvda_flat_biv_to_dft = spqlios_flat_biv_to_dft;
	vt->pvda_flat_biv_idft   = spqlios_flat_biv_idft;

	vt->pvda_svp_apply_dft        = spqlios_svp_apply_dft;
	vt->pvda_svp_apply_dft_to_dft = spqlios_svp_apply_dft_to_dft;

	vt->pvda_new_vec_znx_big    = spqlios_new_vec_znx_big;
	vt->pvda_delete_vec_znx_big = spqlios_delete_vec_znx_big;
	
	vt->pvda_new_biv_vmp_prep     = spqlios_new_biv_vmp_prep;
	vt->pvda_delete_biv_vmp_prep  = spqlios_delete_biv_vmp_prep;
	vt->pvda_flat_biv_to_vmp_prep = spqlios_flat_biv_to_vmp_prep;

	vt->pvda_new_vmp_pmat 			= spqlios_new_vmp_pmat;
	vt->pvda_vmp_prepare_contiguous = spqlios_vmp_prepare_contiguous;
	vt->pvda_delete_vmp_pmat        = spqlios_delete_vmp_pmat;

	vt->pvda_vmp_apply_dft             = spqlios_vmp_apply_dft;
	vt->pvda_vmp_apply_dft_to_dft      = spqlios_vmp_apply_dft_to_dft;
	vt->pvda_vmp_apply_prepared_to_dft = spqlios_vmp_apply_prepared_to_dft;

	vt->pvda_vec_znx_normalize_base2k  = spqlios_vec_znx_normalize_base2k;
	vt->pvda_znx_small_product         = spqlios_znx_small_product;
	vt->pvda_vec_znx_negate            = spqlios_vec_znx_negate;
	vt->pvda_vec_znx_add               = spqlios_vec_znx_add;
	vt->pvda_vec_znx_sub               = spqlios_vec_znx_sub;
	vt->pvda_znx_automorphism          = spqlios_znx_automorphism;
	vt->pvda_vec_znx_automorphism      = spqlios_vec_znx_automorphism;
	vt->pvda_vec_znx_rotate            = spqlios_vec_znx_rotate;
	vt->pvda_module_extract_nn         = spqlios_module_extract_nn;
	vt->pvda_vec_znx_dft_add           = spqlios_vec_znx_dft_add;

	// ! NOTE : These vtable entries are filled elsewhere (rng_openssl.c, cpu_nt.c) :
	// PVDA_RAND_UNIFORM_POW2_F* pvda_rand_uniform_pow2;
	// PVDA_RAND_UNIFORM_POW2_VEC_F* pvda_rand_uniform_pow2_vec;
	// PVDA_RAND_UNIFORM_BINARY_VEC_F* pvda_rand_uniform_binary_vec;
	// PVDA_RAND_UNIFORM_F* pvda_rand_uniform;
	// PVDA_RING_2NTH_ROOT_F* pvda_ring_2nth_root;
}
