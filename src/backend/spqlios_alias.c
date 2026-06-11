#include "spqlios_alias.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#include "rng.h"
#include "utils.h"
#include "vec_znx_arithmetic.h"

MODULE* pvda_new_module_info(uint64_t nn) { return new_module_info(nn, FFT64); }

void pvda_delete_module_info(MODULE* module)
{
	if (!module) return;
	delete_module_info(module);
}

VecUnivDFT* pvda_new_vec_znx_dft(const MODULE* module, uint64_t size) { return (double*)new_vec_znx_dft(module, size); }

void pvda_vec_znx_dft(const MODULE* module, double* res, uint64_t res_size, const PolyBiv* a)
{
	vec_znx_dft(module, (VEC_ZNX_DFT*)res, res_size, a->ptr, a->l, a->stride);
}

void pvda_delete_vec_znx_dft(double* res)
{
	if (!res) return;
	delete_vec_znx_dft((VEC_ZNX_DFT*)res);
}

int64_t* pvda_new_vec_znx_big(const MODULE* module, int64_t size) { return (int64_t*)new_vec_znx_big(module, size); }

int pvda_vec_znx_idft(const MODULE* module, PolyBiv* res, const double* a_dft, uint64_t a_size)
{
	size_t tmp_bytes = vec_znx_idft_tmp_bytes(module);
	uint8_t* tmp     = malloc(tmp_bytes);
	if (tmp_bytes) CHECK_ALLOC(tmp, "ntt120_vec_znx_idft_tmp_bytes_avx failed in vec_znx_idft_p");

	assert(res->stride == res->nn);
	vec_znx_idft(module, (VEC_ZNX_BIG*)res->ptr, res->l, (VEC_ZNX_DFT*)a_dft, a_size, tmp);

cleanup:
	free(tmp);

	return 0;
}

void pvda_delete_vec_znx_big(int64_t* res)
{
	if (!res) return;
	delete_vec_znx_big((VEC_ZNX_BIG*)res);
}

double* pvda_new_svp_ppol(const MODULE* module) { return (PolyUnivDFT*)new_svp_ppol(module); }

void pvda_svp_prepare(const MODULE* module, PolyUnivDFT* ppol, const int64_t* pol)
{
	svp_prepare(module, (SVP_PPOL*)ppol, pol);
}

void pvda_delete_svp_ppol(double* pmat)
{
	if (!pmat) return;
	delete_svp_ppol(((SVP_PPOL*)pmat));
}

void pvda_svp_apply_dft(const MODULE* module, const double* res, uint64_t res_size, const PolyUnivDFT* prepared_pol,
                        const PolyBiv* a)
{
	svp_apply_dft(module, (VEC_ZNX_DFT*)res, res_size, (SVP_PPOL*)prepared_pol, a->ptr, a->l, a->stride);
}

double* pvda_new_vmp_pmat(const MODULE* module, uint64_t nrows, uint64_t ncols)
{
	return (double*)new_vmp_pmat(module, nrows, ncols);
}

int pvda_vmp_prepare_contiguous(const MODULE* module, double* pmat, const int64_t* mat, uint64_t nrows, uint64_t ncols)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vmp_prepare_contiguous_tmp_bytes(module, nrows, ncols));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vmp_prepare_contiguous_p");

	vmp_prepare_contiguous(module, (VMP_PMAT*)pmat, mat, nrows, ncols, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

void pvda_delete_vmp_pmat(double* pmat)
{
	if (!pmat) return;
	delete_vmp_pmat(((VMP_PMAT*)pmat));
}

int pvda_vmp_apply_dft(const MODULE* module, double* res, uint64_t res_size, const PolyBiv* a, const MatBivDFT* pmat,
                       uint64_t nrows, uint64_t ncols)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vmp_apply_dft_tmp_bytes(module, res_size, a->l, nrows, ncols));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vmp_apply_dft_p");

	vmp_apply_dft(module, (VEC_ZNX_DFT*)res, res_size, a->ptr, a->l, a->stride, (VMP_PMAT*)pmat, nrows, ncols,
	              tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

int pvda_vmp_apply_dft_to_dft(const MODULE* module, VecBivDFT* res, const uint64_t res_size, const VecBivDFT* a_dft,
                              uint64_t a_size, const MatBivDFT* pmat, const uint64_t nrows, const uint64_t ncols)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vmp_apply_dft_to_dft_tmp_bytes(module, res_size, a_size, nrows, ncols));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vmp_apply_dft_to_dft_p");

	vmp_apply_dft_to_dft(module, (VEC_ZNX_DFT*)res, res_size, (VEC_ZNX_DFT*)a_dft, a_size, (VMP_PMAT*)pmat, nrows,
	                     ncols, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

int pvda_vmp_apply_prepared_to_dft(const MODULE* module, VecBivDFT* res, const uint64_t res_size,
                                   const VecBivDFT* a_dft, uint64_t a_size, const MatBivDFT* pmat, const uint64_t nrows,
                                   const uint64_t ncols)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vmp_apply_prepared_to_dft_tmp_bytes(module, res_size, a_size, nrows, ncols));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vmp_apply_dft_to_dft_p");

	vmp_apply_prepared_to_dft(module, (VEC_ZNX_DFT*)res, res_size, (VEC_ZNX_DFT*)a_dft, a_size, (VMP_PMAT*)pmat, nrows,
	                          ncols, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

int pvda_vec_znx_normalize_base2k(const MODULE* module, uint64_t log2_base2k, PolyBiv* res, const PolyBiv* a)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vec_znx_normalize_base2k_tmp_bytes(module));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vec_znx_normalize_base2k_p");

	vec_znx_normalize_base2k(module, log2_base2k, res->ptr, res->l, res->stride, a->ptr, a->l, a->stride, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

int pvda_znx_small_product(const MODULE* module, PolyUniv* res, const PolyUniv* a, const PolyUniv* b)
{
	int status         = -1;
	size_t tmp_size    = znx_small_single_product_tmp_bytes(module);
	uint8_t* tmp_space = malloc(tmp_size);
	if (tmp_size) CHECK_ALLOC(tmp_space, "failed temp space alloc for polynomial multiplication");
	znx_small_single_product(module, res, a, b, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}

int pvda_vec_znx_negate(const MODULE* module, PolyBiv* res, const PolyBiv* a)
{
	vec_znx_negate(module, res->ptr, res->l, res->stride, a->ptr, a->l, a->stride);
	return 0;
}

uint64_t pvda_module_extract_nn(const MODULE* module) { return module_get_n(module); }

int pvda_vec_znx_add(const MODULE* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b)
{
	vec_znx_add(module, res->ptr, res->l, res->stride, a->ptr, a->l, a->stride, b->ptr, b->l, b->stride);
	return 0;
}
int pvda_vec_znx_sub(const MODULE* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b)
{
	vec_znx_sub(module, res->ptr, res->l, res->stride, a->ptr, a->l, a->stride, b->ptr, b->l, b->stride);
	return 0;
}

int pvda_znx_automorphism(const MODULE* module, const int64_t p, PolyUniv* res, const PolyUniv* a)
{
	vec_znx_automorphism(module, p, res, 1, 0, a, 1, 0);
	return 0;
}
int pvda_vec_znx_automorphism(const MODULE* module, const int64_t p, PolyBiv* res, const PolyBiv* a)
{
	vec_znx_automorphism(module, p, res->ptr, res->l, res->stride, a->ptr, a->l, a->stride);
	return 0;
}

int pvda_vec_znx_rotate(const MODULE* module, const int64_t p, PolyBiv* res, const PolyBiv* a)
{
	vec_znx_rotate(module, p, res->ptr, res->l, res->stride, a->ptr, a->l, a->stride);
	return 0;
}

int pvda_vmp_prepare_vec(const MODULE* module, double* pvec, uint64_t nrows, const PolyBiv* a)
{
	int status = -1;

	void* tmp_space = aligned_alloc(64, cnv_prepare_right_contiguous_tmp_bytes(module, nrows, a->l));
	CHECK_ALLOC(tmp_space, "tmp space alloc failed in pvda_convolution_prepare");

	cnv_prepare_right_contiguous(module, (CNV_PVEC_R*)pvec, nrows, a->ptr, a->l, a->stride, tmp_space);
	status = 0;
cleanup:
	free(tmp_space);
	return status;
}
/*
int pvda_vec_rnx_negate(const MODULE* module, double* res, uint64_t res_size, uint64_t res_sl, const double* a,
                        uint64_t a_size, uint64_t a_sl)
{
    vec_rnx_negate(module, res, res_size, res_sl, a, a_size, a_sl);
    return 1;
}*/
