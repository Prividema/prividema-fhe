#ifndef PVDA_BACKEND_ARITH
#define PVDA_BACKEND_ARITH

#include "backend.h"
#include "pvda_types.h"

/**
 * @file backend_arithmetic.h
 *
 * @brief Arithmetic functions of the adapter.
 *
 * This file contains all the arithmetic functions from the backend we need on the library side.
 *
 * Every function takes library-side data (PolyUnivZnx, PolyBiv) or opaque `BE_*` handles. The adapter converts
 * library-side data into the input format of the selected backend (for spqlios, a flat strided int64 array) and
 * casts the `BE_*` handles to the true backend types. No backend-specific type leaks through this interface.
 *
 * Univariate PolyUnivZnx arguments are \N contiguous coefficients, \N being the ring dimension of the module.
 *
 * Bivariate PolyBiv arguments are read as a strided vector of nb_limbs \ZnX polynomials (see PolyBiv). Unless stated
 * otherwise, their nn must match the ring dimension of the module.
 *
 * All functions returning an int return 0 on success and -1 on failure.
 */

// =============================================
// |                                           |
// |              BE_PolyUnivZnxDFT            |
// |                                           |
// =============================================

BE_PolyUnivZnxDFT* pvda_new_univ_znx_dft(const PvdaBackend* module);
void pvda_delete_univ_znx_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* pol);

/** @brief DFT of a univariate polynomial. */
int pvda_univ_znx_to_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* res, const PolyUnivZnx* a);

/**
 * @brief Inverse DFT of a univariate polynomial into coefficient form.
 *
 * Output coefficients are NOT normalized (no carry propagation): each coefficient is the exact
 * int64 value of the inverse transform. Callers that need base-2^k limbs must normalize.
 * Any backend-specific "big" intermediate layout is handled inside the backend.
 */
int pvda_univ_znx_idft(const PvdaBackend* module, PolyUnivZnx* res, const BE_PolyUnivZnxDFT* pol);

/**
 * @brief View of a univariate DFT polynomial as a 1-limb bivariate DFT polynomial. Does not allocate or copy.
 *
 * The returned pointer aliases pol: it is valid as long as pol is and must not be deleted.
 */
BE_PolyBivDFT* pvda_univ_dft_as_biv_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* pol);

/** @brief Const version of pvda_univ_dft_as_biv_dft. */
const BE_PolyBivDFT* pvda_univ_dft_as_biv_dft_const(const PvdaBackend* module, const BE_PolyUnivZnxDFT* pol);

// =============================================
// |                                           |
// |             BE_PolyUnivSvpPrep            |
// |                                           |
// =============================================

BE_PolyUnivSvpPrep* pvda_new_univ_svp_prep(const PvdaBackend* module);
void pvda_delete_univ_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* pol);

/** @brief Prepares a univariate polynomial for scalar-vector products. */
int pvda_univ_znx_to_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* res, const PolyUnivZnx* a);

// =============================================
// |                                           |
// |                BE_PolyBivDFT              |
// |                                           |
// =============================================

BE_PolyBivDFT* pvda_new_biv_dft(const PvdaBackend* module, uint64_t size);
void pvda_delete_biv_dft(const PvdaBackend* module, BE_PolyBivDFT* res);

/**
 * @brief DFT of the a->nb_limbs polynomials of a into the res_size DFT polynomials of res.
 *
 * If res_size > a->nb_limbs, the extra DFT polynomials are set to zero.
 */
int pvda_biv_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size, const PolyBiv* a);

/**
 * @brief Inverse DFT of a vector of pol_size DFT polynomials into res->nb_limbs coefficient-form polynomials.
 *
 * Output limbs are NOT normalized (no carry propagation): each coefficient is the exact
 * int64 value of the inverse transform. Callers that need base-2^k limbs must normalize.
 * Any backend-specific "big" intermediate layout is handled inside the backend.
 *
 * @pre res->stride == res->nn (returns -1 otherwise)
 */
int pvda_biv_idft(const PvdaBackend* module, PolyBiv* res, const BE_PolyBivDFT* pol, uint64_t pol_size);

/**
 * @brief View of a DFT vector starting at its offset'th DFT polynomial. Does not allocate or copy.
 *
 * The returned pointer aliases vec: it is valid as long as vec is, must not be deleted, and holds the
 * DFT polynomials offset, offset+1, ... of vec (as many as vec has after offset).
 *
 * The view is contiguous: DFT vectors have no stride. In a GLWE DFT, limbs are stored in the same interleaved
 * order as in coefficient domain (limb i of element j at index i * (k+1) + j), so a view at offset i * (k+1)
 * gives limb i of all k+1 elements, and a 1-limb view gives one limb of one element. An element of the GLWE
 * (all limbs of a_j) is strided and cannot be viewed this way.
 *
 * @param module The computation backend
 * @param vec    The DFT vector to view
 * @param offset Index, in DFT polynomials, of the first polynomial of the view
 *
 * @return A pointer to the offset'th DFT polynomial of vec
 */
BE_PolyBivDFT* pvda_biv_dft_at(const PvdaBackend* module, BE_PolyBivDFT* vec, uint64_t offset);

/** @brief Const version of pvda_biv_dft_at. */
const BE_PolyBivDFT* pvda_biv_dft_at_const(const PvdaBackend* module, const BE_PolyBivDFT* vec, uint64_t offset);

// =============================================
// |                                           |
// |                    SVP                    |
// |                                           |
// =============================================

/** @brief res = DFT(prepared_pol * a), polynomial by polynomial, over res_size DFT polynomials. */
int pvda_svp_apply_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size,
                       const BE_PolyUnivSvpPrep* prepared_pol, const PolyBiv* a);

/** @brief res = ppol * a in the DFT domain, polynomial by polynomial. */
int pvda_svp_apply_dft_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size,
                              const BE_PolyUnivSvpPrep* ppol, const BE_PolyBivDFT* a, uint64_t a_size);

// =============================================
// |                                           |
// |             BE_PolyBivVmpPrep             |
// |                                           |
// =============================================

BE_PolyBivVmpPrep* pvda_new_biv_vmp_prep(const PvdaBackend* module, uint64_t biv_l);
void pvda_delete_biv_vmp_prep(const PvdaBackend* module, BE_PolyBivVmpPrep* biv_prep);

/** @brief Prepares the polynomials of a as a vector (of res_size rows) for vector-matrix products. */
int pvda_biv_to_vmp_prep(const PvdaBackend* module, BE_PolyBivVmpPrep* res, uint64_t res_size, const PolyBiv* a);

// =============================================
// |                                           |
// |              BE_PolyMatVmpPrep            |
// |                                           |
// =============================================

BE_PolyMatVmpPrep* pvda_new_vmp_pmat(const PvdaBackend* module, uint64_t nrows, uint64_t ncols);
void pvda_delete_vmp_pmat(const PvdaBackend* module, BE_PolyMatVmpPrep* pmat);

/**
 * @brief Prepares a matrix of nrows * ncols polynomials for vector-matrix products.
 *
 * The matrix is given row-major as the polynomials of mat: polynomial (i, j) is limb i * ncols + j of mat.
 *
 * @pre mat->nb_limbs == nrows * ncols (returns -1 otherwise)
 */
int pvda_biv_to_vmp_pmat(const PvdaBackend* module, BE_PolyMatVmpPrep* pmat, const PolyBiv* mat, uint64_t nrows,
                         uint64_t ncols);

// =============================================
// |                                           |
// |                    VMP                    |
// |                                           |
// =============================================

int pvda_vmp_apply_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size, const PolyBiv* a,
                       const BE_PolyMatVmpPrep* pmat, uint64_t nrows, uint64_t ncols);
int pvda_vmp_apply_dft_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, const uint64_t res_size,
                              const BE_PolyBivDFT* a_dft, uint64_t a_size, const BE_PolyMatVmpPrep* pmat,
                              const uint64_t nrows, const uint64_t ncols);
int pvda_vmp_apply_prepared_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, const uint64_t res_size,
                                   const BE_PolyBivVmpPrep* a_prep, uint64_t a_size, const BE_PolyMatVmpPrep* pmat,
                                   const uint64_t nrows, const uint64_t ncols);

// =============================================
// |                                           |
// |                 Arithmetic                |
// |                                           |
// =============================================

/** @brief res = a with carry propagated across limbs in base 2^log2_base2k. */
int pvda_biv_normalize_base2k(const PvdaBackend* module, uint64_t log2_base2k, PolyBiv* res, const PolyBiv* a);

/** @brief res = a * b in \ZnX (exact for small coefficients). */
int pvda_univ_znx_small_product(const PvdaBackend* module, PolyUnivZnx* res, const PolyUnivZnx* a,
                                const PolyUnivZnx* b);

/** @brief res = -a, limb by limb. */
int pvda_biv_negate(const PvdaBackend* module, PolyBiv* res, const PolyBiv* a);

/** @brief res = a + b, limb by limb. Missing limbs of the shortest input are read as zero. */
int pvda_biv_add(const PvdaBackend* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b);

/** @brief res = a - b, limb by limb. Missing limbs of the shortest input are read as zero. */
int pvda_biv_sub(const PvdaBackend* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b);

/** @brief res = a(X^p) for a univariate polynomial. */
int pvda_univ_znx_automorphism(const PvdaBackend* module, const int64_t p, PolyUnivZnx* res, const PolyUnivZnx* a);

/** @brief res = a(X^p), limb by limb. */
int pvda_biv_automorphism(const PvdaBackend* module, const int64_t p, PolyBiv* res, const PolyBiv* a);

/** @brief res = X^p * a, limb by limb. */
int pvda_biv_rotate(const PvdaBackend* module, const int64_t p, PolyBiv* res, const PolyBiv* a);

/** @brief res = a + b in the DFT domain. Missing DFT polynomials of the shortest input are read as zero. */
int pvda_biv_dft_add(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_n_limbs, const BE_PolyBivDFT* a,
                     uint64_t a_n_limbs, const BE_PolyBivDFT* b, uint64_t b_n_limbs);
#endif
