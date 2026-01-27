#include "spqlios_alias.h"

#include "vec_znx_arithmetic_private.h"

MODULE* new_module_info_p(uint64_t N) {
  return new_module_info(N, FFT64);
}

void delete_module_info_p(MODULE* module){
  delete_module_info(module);
}

double* new_vec_znx_dft_p(const MODULE* module, int64_t size) {
  return (double*)new_vec_znx_dft(module, size);
}

void vec_znx_dft_p(const MODULE* module, double* res, int64_t res_size,
                   const int64_t* a, int64_t a_size, int64_t a_sl) {
  vec_znx_dft(module, (VEC_ZNX_DFT*)res, res_size, a, a_size, a_sl);
}

void delete_vec_znx_dft_p(double* res) {
  delete_vec_znx_dft((VEC_ZNX_DFT*)res);
}

int64_t* new_vec_znx_big_p(const MODULE* module, int64_t size) {
  return (int64_t*)new_vec_znx_big(module, size);
}

void vec_znx_idft_p(const MODULE* module, int64_t* res, int64_t res_size,
                    const double* a_dft, int64_t a_size) {
  uint8_t* tmp = NULL;
  vec_znx_idft(module, (VEC_ZNX_BIG*)res, res_size, (VEC_ZNX_DFT*)a_dft, a_size,
               tmp);
}

void delete_vec_znx_big_p(int64_t* res) {
  delete_vec_znx_big((VEC_ZNX_BIG*)res);
}

double* new_svp_ppol_p(MODULE* module) {
  return (PolyUnivDFT*)new_svp_ppol(module);
}

void svp_prepare_p(const MODULE* module, PolyUnivDFT* ppol,
                   const int64_t* pol) {
  svp_prepare(module, (SVP_PPOL*)ppol, pol);
}

void delete_svp_ppol_p(double* pmat) { delete_svp_ppol(((SVP_PPOL*)pmat)); }

void svp_apply_dft_p(const MODULE* module, const double* res, int64_t res_size,
                     const PolyUnivDFT* ppol, const int64_t* a, int64_t a_size,
                     int64_t a_sl) {
    svp_apply_dft(module, (VEC_ZNX_DFT*)res, res_size, (SVP_PPOL*)ppol, a, a_size, a_sl);
}

double* new_vmp_pmat_p(const MODULE* module, uint64_t nrows, uint64_t ncols) {
    return (double*)new_vmp_pmat(module, nrows, ncols);
}

void vmp_prepare_contiguous_p(const MODULE* module, double* pmat,
                              const int64_t* mat, uint64_t nrows,
                              uint64_t ncols)
{
    uint8_t* tmp_space = malloc(vmp_prepare_contiguous_tmp_bytes(module, nrows, ncols));
    vmp_prepare_contiguous(module, (VMP_PMAT*)pmat, mat, nrows, ncols, tmp_space);
    free(tmp_space);
}

void delete_vmp_pmat_p(double* pmat) { delete_vmp_pmat(((VMP_PMAT*)pmat)); }

void vmp_apply_dft_p(const MODULE* module, double* res, int64_t res_size,
                     const int64_t* a, int64_t a_size, int64_t a_sl,
                     const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols)
{
    uint8_t* tmp_space = malloc(vmp_apply_dft_tmp_bytes(module, res_size, a_size, nrows, ncols));
    vmp_apply_dft(module, (VEC_ZNX_DFT*)res, res_size, a, a_size, a_sl,
                  (VMP_PMAT*)pmat, nrows, ncols, tmp_space);
    free(tmp_space);
}

void vmp_apply_dft_to_dft_p(const MODULE* module, VecBivDFT* res,
                            const uint64_t res_size, const VecBivDFT* a_dft,
                            uint64_t a_size, const MatBivDFT* pmat,
                            const uint64_t nrows, const uint64_t ncols) {
    uint8_t* tmp_space = malloc(vmp_apply_dft_to_dft_tmp_bytes(module, res_size, a_size, nrows, ncols));
    vmp_apply_dft_to_dft(module, (VEC_ZNX_DFT*)res, res_size, (VEC_ZNX_DFT*)a_dft,
                         a_size, (VMP_PMAT*)pmat, nrows, ncols, tmp_space);
    free(tmp_space);
}

void vec_znx_normalize_base2k_p(const MODULE* module, uint64_t log2_base2k,
                                int64_t* res, int64_t res_size,
                                int64_t res_sl, const int64_t* a,
                                int64_t a_size, int64_t a_sl) {
    uint8_t* tmp = malloc(vec_znx_normalize_base2k_tmp_bytes(module));
    vec_znx_normalize_base2k(module, log2_base2k, res, res_size, res_sl, a,
                             a_size, a_sl, tmp);
    free(tmp);
}