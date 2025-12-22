#include "structure_alias.h"


//! DEFINE SPQLIOS ALIAS (begin) 

/** @brief opaque structure that describr the modules (ZnX,TnX) and the hardware */
typedef struct module_info_t MODULE;
/** @brief opaque type that represents a prepared matrix */
typedef struct vmp_pmat_t VMP_PMAT;
/** @brief opaque type that represents a vector of znx in DFT space */
typedef struct vec_znx_dft_t VEC_ZNX_DFT;
/** @brief opaque type that represents a vector of znx in large coeffs space */
typedef struct vec_znx_bigcoeff_t VEC_ZNX_BIG;
/** @brief opaque type that represents a prepared scalar vector product */
typedef struct svp_ppol_t SVP_PPOL;
/** @brief opaque type that represents a prepared left convolution vector product */
typedef struct cnv_pvec_l_t CNV_PVEC_L;
/** @brief opaque type that represents a prepared right convolution vector product */
typedef struct cnv_pvec_r_t CNV_PVEC_R;

//! The p in "*_p" stands for PRIVIDEMA (begin)

double* new_vec_znx_dft_p(const MODULE* module,  // N
                          uint64_t size
);

void vec_znx_dft_p(const MODULE* module,                             // N
                   double* res, uint64_t res_size,                   // res
                   const int64_t* a, uint64_t a_size, uint64_t a_sl  // a
);

void delete_vec_znx_dft_p(double* res);

int64_t* new_vec_znx_big_p(const MODULE* module,  // N
                           uint64_t size
);

void vec_znx_idft_p(const MODULE* module,                  // N
                    int64_t* res, uint64_t res_size,       // res
                    const double* a_dft, uint64_t a_size,  // a
                    uint8_t* tmp                           // useless
);

void delete_vec_znx_big_p(int64_t* res);

double* new_svp_ppol_p(MODULE* module);

void svp_prepare_p(const MODULE* module,  // N
                        PolyUnivDFT* ppol,        // output
                        const int64_t* pol     // a
);

void svp_apply_dft_p(const MODULE* module,                             // N
                     const double* res, uint64_t res_size,             // output
                     const PolyUnivDFT* ppol,                          // prepared pol
                     const int64_t* a, uint64_t a_size, uint64_t a_sl  // a
);

void delete_svp_ppol_p(double* res);

double* new_vmp_pmat_p(const MODULE* module, uint64_t nrows, uint64_t ncols);

void vmp_prepare_contiguous_p(const MODULE* module,                                // N
                                   double* pmat,                                      // output
                                   const int64_t* mat, uint64_t nrows, uint64_t ncols,  // a
                                   uint8_t* tmp_space                                   // scratch space
);

void vmp_apply_dft_p(const MODULE* module,                                   // N
                     double* res, uint64_t res_size,                         // res
                     const int64_t* a, uint64_t a_size, uint64_t a_sl,       // a
                     const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols,  // prep matrix
                     uint8_t* tmp_space                                      // scratch space
);

void vmp_apply_dft_to_dft_p(const MODULE* module,                       // N
                          double* res, const uint64_t res_size,  // res
                          const double* a_dft, uint64_t a_size,  // a
                          const MatBivDFT* pmat, const uint64_t nrows,
                          const uint64_t ncols,  // prep matrix
                          uint8_t* tmp_space     // scratch space (a_size*sizeof(reim4) bytes)
);

void delete_vmp_pmat_p(double* pmat);