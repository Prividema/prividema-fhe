#include "maths_structures.h"
#include "vec_znx_arithmetic_private.h"

// =============================================
// |                                           |
// |      Aliases for spqlios structures       |
// |                                           |
// =============================================

typedef enum module_type_t MODULE_TYPE;

/// @brief Opaque structure that describe the modules (\f$\mathbb{Z}_n[X]\f$, \f$\mathbb{T}_n[X]\f$) and the hardware.
typedef struct module_info_t MODULE;

/// @brief Opaque type that represents a prepared matrix.
typedef struct vmp_pmat_t VMP_PMAT;

/// @brief Opaque type that represents a vector of \f$\mathbb{Z}_n[X]\f$ in the DFT domain.
typedef struct vec_znx_dft_t VEC_ZNX_DFT;

/// @brief Opaque type that represents a vector of \f$\mathbb{Z}_n[X]\f$ in large coeffs space.
typedef struct vec_znx_bigcoeff_t VEC_ZNX_BIG;

/// @brief Opaque type that represents a prepared scalar vector product.
typedef struct svp_ppol_t SVP_PPOL;

/// @brief Opaque type that represents a prepared left convolution vector product.
typedef struct cnv_pvec_l_t CNV_PVEC_L;

/// @brief Opaque type that represents a prepared right convolution vector product.
typedef struct cnv_pvec_r_t CNV_PVEC_R;

// =============================================
// |                                           |
// |             spqlios Functions             |
// |                                           |
// |       The following function are all      |
// |      named f_p where f is the name of     |
// |           a function in sqplios.          |
// |                                           |
// |  We additionnally allocate the memory of  |
// |                the structures.            |
// |                                           |
// =============================================

MODULE* new_module_info_p(uint64_t N);

void delete_module_info_p(MODULE* module);

double* new_vec_znx_dft_p(const MODULE* module, int64_t size);

void vec_znx_dft_p(const MODULE* module,
                   double* res, int64_t res_size,
                   const int64_t* a, int64_t a_size, int64_t a_sl);

void delete_vec_znx_dft_p(double* res);

int64_t* new_vec_znx_big_p(const MODULE* module, int64_t size);

int vec_znx_idft_p(const MODULE* module,
                   int64_t* res, int64_t res_size,
                   const double* a_dft, int64_t a_size);

void delete_vec_znx_big_p(int64_t* res);

double* new_svp_ppol_p(const MODULE* module);

void svp_prepare_p(const MODULE* module,
                   PolyUnivDFT* prepared_pol,
                   const int64_t* pol);

void svp_apply_dft_p(const MODULE* module,
                     const double* res, int64_t res_size,
                     const PolyUnivDFT* prepared_pol,
                     const int64_t* a, int64_t a_size, int64_t a_sl);

void delete_svp_ppol_p(double* res);

double* new_vmp_pmat_p(const MODULE* module, uint64_t nrows, uint64_t ncols);

int vmp_prepare_contiguous_p(const MODULE* module,
                              double* pmat,
                              const int64_t* mat, uint64_t nrows, uint64_t ncols);

int vmp_apply_dft_p(const MODULE* module,
                     double* res, int64_t res_size,
                     const int64_t* a, int64_t a_size, int64_t a_sl,
                     const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols);

int vmp_apply_dft_to_dft_p(const MODULE* module,
                           VecBivDFT* res, const uint64_t res_size,
                           const VecBivDFT* a_dft, uint64_t a_size,
                           const MatBivDFT* pmat, const uint64_t nrows, const uint64_t ncols);

void delete_vmp_pmat_p(double* pmat);

int vec_znx_normalize_base2k_p(const MODULE* module,
                                uint64_t log2_base2k,
                                int64_t* res, int64_t res_size, int64_t res_sl,
                                const int64_t* a, int64_t a_size, int64_t a_sl);