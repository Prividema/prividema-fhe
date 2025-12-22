#include "spqlios_alias.h"
#include "vec_znx_arithmetic.h"

double* new_vec_znx_dft_p(const MODULE* module,  // N
                          uint64_t size
){
    return (double*)new_vec_znx_dft(module, size);
}

void vec_znx_dft_p(const MODULE* module,                             // N
                   double* res, uint64_t res_size,                   // res
                   const int64_t* a, uint64_t a_size, uint64_t a_sl  // a
){
    vec_znx_dft(module, (VEC_ZNX_DFT*)res, res_size, a, a_size, a_sl);
}

void delete_vec_znx_dft_p(double* res){
    delete_vec_znx_dft((VEC_ZNX_DFT*)res);
}

int64_t* new_vec_znx_big_p(const MODULE* module,  // N
                           uint64_t size
){
    return (int64_t*)new_vec_znx_big(module, size);
}

void vec_znx_idft_p(const MODULE* module,                  // N
                    int64_t* res, uint64_t res_size,       // res
                    const double* a_dft, uint64_t a_size,  // a
                    uint8_t* tmp                           // useless
){
    vec_znx_idft(module, (VEC_ZNX_BIG*)res, res_size, (VEC_ZNX_DFT*)a_dft, a_size, tmp);
}

void delete_vec_znx_big_p(int64_t* res){
    delete_vec_znx_big((VEC_ZNX_BIG*)res);
}

double* new_svp_ppol_p(MODULE* module){
    return (PolyUnivDFT*)new_svp_ppol(module);
}

void svp_prepare_p(const MODULE* module,  // N
                        PolyUnivDFT* ppol,        // output
                        const int64_t* pol     // a
){
    svp_prepare(module, (SVP_PPOL*)ppol, pol);
}

void delete_svp_ppol_p(double* pmat){
    delete_svp_ppol(((SVP_PPOL*)pmat));
}

void svp_apply_dft_p(const MODULE* module,                             // N
                     const double* res, uint64_t res_size,             // output
                     const PolyUnivDFT* ppol,                          // prepared pol
                     const int64_t* a, uint64_t a_size, uint64_t a_sl  // a
){
    svp_apply_dft(module, (VEC_ZNX_DFT*)res, res_size, (SVP_PPOL*)ppol, a, a_size, a_sl);
}

double* new_vmp_pmat_p(const MODULE* module,           // N
                       uint64_t nrows, uint64_t ncols  // dimensions
){
    return (double*)new_vmp_pmat(module, nrows, ncols);
}

void vmp_prepare_contiguous_p(const MODULE* module,                                // N
                              double* pmat,                                      // output
                              const int64_t* mat, uint64_t nrows, uint64_t ncols,  // a
                              uint8_t* tmp_space                                   // scratch space
){
    vmp_prepare_contiguous(module, (VMP_PMAT*)pmat, mat, nrows, ncols, tmp_space);
}

void delete_vmp_pmat_p(double* pmat){
    delete_vmp_pmat(((VMP_PMAT*)pmat));
}

void vmp_apply_dft_p(const MODULE* module,                                   // N
                     double* res, uint64_t res_size,                         // res
                     const int64_t* a, uint64_t a_size, uint64_t a_sl,       // a
                     const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols,  // prep matrix
                     uint8_t* tmp_space                                      // scratch space
){
    vmp_apply_dft(module, (VEC_ZNX_DFT*)res, res_size, a, a_size, a_sl, (VMP_PMAT*)pmat,nrows, ncols, tmp_space);
}

void vmp_apply_dft_to_dft_p(const MODULE* module,                       // N
                          double* res, const uint64_t res_size,  // res
                          const double* a_dft, uint64_t a_size,  // a
                          const MatBivDFT* pmat, const uint64_t nrows,
                          const uint64_t ncols,  // prep matrix
                          uint8_t* tmp_space     // scratch space (a_size*sizeof(reim4) bytes)
){
    vmp_apply_dft_to_dft(module, (VEC_ZNX_DFT*)res, res_size, (VEC_ZNX_DFT*)a_dft, a_size, (VMP_PMAT*)pmat, nrows, ncols, tmp_space);
}