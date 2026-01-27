#include "structure_alias.h"


//! DEFINE SPQLIOS ALIAS (begin) 

typedef enum module_type_t MODULE_TYPE;
/** @brief opaque structure that describe the modules (Zn[X],TnX) and the hardware */
typedef struct module_info_t MODULE;
/** @brief opaque type that represents a prepared matrix */
typedef struct vmp_pmat_t VMP_PMAT;
/** @brief opaque type that represents a vector of Zn[X] in DFT space */
typedef struct vec_znx_dft_t VEC_ZNX_DFT;
/** @brief opaque type that represents a vector of Zn[X] in large coeffs space */
typedef struct vec_znx_bigcoeff_t VEC_ZNX_BIG;
/** @brief opaque type that represents a prepared scalar vector product */
typedef struct svp_ppol_t SVP_PPOL;
/** @brief opaque type that represents a prepared left convolution vector product */
typedef struct cnv_pvec_l_t CNV_PVEC_L;
/** @brief opaque type that represents a prepared right convolution vector product */
typedef struct cnv_pvec_r_t CNV_PVEC_R;

//! The p in "*_p" stands for PRIVIDEMA (begin)

/**
 * @brief Creates a new module holding N referring to Z_N[X] and we suppose that we use the FFT64 Discrete Fourier Transform method.
 * 
 * @param N 
 * @return MODULE* 
 */
MODULE* new_module_info_p(uint64_t N
);

void delete_module_info_p(MODULE* module);

double* new_vec_znx_dft_p(const MODULE* module,  // N
                          int64_t size
);

void vec_znx_dft_p(const MODULE* module,                             // N
                   double* res, int64_t res_size,                   // res
                   const int64_t* a, int64_t a_size, int64_t a_sl  // a
);

void delete_vec_znx_dft_p(double* res);

int64_t* new_vec_znx_big_p(const MODULE* module,  // N
                           int64_t size
);

void vec_znx_idft_p(const MODULE* module,                  // N
                    int64_t* res, int64_t res_size,       // res
                    const double* a_dft, int64_t a_size   // a
);

void delete_vec_znx_big_p(int64_t* res);

double* new_svp_ppol_p(MODULE* module);

void svp_prepare_p(const MODULE* module,  // N
                   PolyUnivDFT* ppol,     // output
                   const int64_t* pol     // a
);

void svp_apply_dft_p(const MODULE* module,                             // N
                     const double* res, int64_t res_size,             // output
                     const PolyUnivDFT* ppol,                          // prepared pol
                     const int64_t* a, int64_t a_size, int64_t a_sl  // a
);

void delete_svp_ppol_p(double* res);

double* new_vmp_pmat_p(const MODULE* module, uint64_t nrows, uint64_t ncols);

void vmp_prepare_contiguous_p(const MODULE* module,                               // N
                              double* pmat,                                       // output
                              const int64_t* mat, uint64_t nrows, uint64_t ncols  // a
);

void vmp_apply_dft_p(const MODULE* module,                                   // N
                     double* res, int64_t res_size,                         // res
                     const int64_t* a, int64_t a_size, int64_t a_sl,       // a
                     const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols   // prep matrix
);

void vmp_apply_dft_to_dft_p(const MODULE* module,                       // N
                            double* res, const int64_t res_size,       // res
                            const double* a_dft, int64_t a_size,       // a
                            const MatBivDFT* pmat, const uint64_t nrows,
                            const uint64_t ncols                        // prep matrix
);

void delete_vmp_pmat_p(double* pmat);

void vec_znx_normalize_base2k_p(const MODULE* module,                              // N
                                uint64_t log2_base2k,                              // output base 2^K
                                int64_t* res, int64_t res_size, int64_t res_sl,  // res
                                const int64_t* a, int64_t a_size, int64_t a_sl   //a
);