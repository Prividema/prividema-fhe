#include "spqlios_alias.h"

#include "utils.h"

MODULE* pvda_new_module_info(uint64_t N) { return new_module_info(N, FFT64); }

void pvda_delete_module_info(MODULE* module)
{
	if (!module) return;
	delete_module_info(module);
}

double* pvda_new_vec_znx_dft(const MODULE* module, int64_t size) { return (double*)new_vec_znx_dft(module, size); }

void pvda_vec_znx_dft(const MODULE* module, double* res, int64_t res_size, const int64_t* a, int64_t a_size,
                      int64_t a_sl)
{
	vec_znx_dft(module, (VEC_ZNX_DFT*)res, res_size, a, a_size, a_sl);
}

void pvda_delete_vec_znx_dft(double* res)
{
	if (!res) return;
	delete_vec_znx_dft((VEC_ZNX_DFT*)res);
}

int64_t* pvda_new_vec_znx_big(const MODULE* module, int64_t size) { return (int64_t*)new_vec_znx_big(module, size); }

int pvda_vec_znx_idft(const MODULE* module, int64_t* res, int64_t res_size, const double* a_dft, int64_t a_size)
{
	uint8_t* tmp = NULL;
	if (module->module_type == NTT120)
	{
		tmp = malloc(ntt120_vec_znx_idft_tmp_bytes_avx(module));
		CHECK_ALLOC(tmp, "ntt120_vec_znx_idft_tmp_bytes_avx failed in vec_znx_idft_p");
	}

	vec_znx_idft(module, (VEC_ZNX_BIG*)res, res_size, (VEC_ZNX_DFT*)a_dft, a_size, tmp);

cleanup:
	if (module->module_type == NTT120) free(tmp);

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

void pvda_svp_apply_dft(const MODULE* module, const double* res, int64_t res_size, const PolyUnivDFT* ppol,
                        const int64_t* a, int64_t a_size, int64_t a_sl)
{
	svp_apply_dft(module, (VEC_ZNX_DFT*)res, res_size, (SVP_PPOL*)ppol, a, a_size, a_sl);
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

int pvda_vmp_apply_dft(const MODULE* module, double* res, int64_t res_size, const int64_t* a, int64_t a_size,
                       int64_t a_sl, const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vmp_apply_dft_tmp_bytes(module, res_size, a_size, nrows, ncols));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vmp_apply_dft_p");

	vmp_apply_dft(module, (VEC_ZNX_DFT*)res, res_size, a, a_size, a_sl, (VMP_PMAT*)pmat, nrows, ncols, tmp_space);

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

int pvda_vec_znx_normalize_base2k(const MODULE* module, uint64_t log2_base2k, int64_t* res, int64_t res_size,
                                  int64_t res_sl, const int64_t* a, int64_t a_size, int64_t a_sl)
{
	int status         = -1;
	uint8_t* tmp_space = malloc(vec_znx_normalize_base2k_tmp_bytes(module));
	CHECK_ALLOC(tmp_space, "tmp_space's malloc failed in vec_znx_normalize_base2k_p");

	vec_znx_normalize_base2k(module, log2_base2k, res, res_size, res_sl, a, a_size, a_sl, tmp_space);

	status = 0;
cleanup:
	free(tmp_space);
	return status;
}
