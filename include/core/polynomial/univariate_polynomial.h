#ifndef UNIVARIATE_POLYNOMIAL_H
#define UNIVARIATE_POLYNOMIAL_H

#include <stdint.h>

#include "backend.h"
#include "backend_types.h"
#include "glwe_params.h"

/**
 * @file univariate_polynomial.h
 *
 * @brief Operations on univariate polynomials.
 * 
 * This file contains everything on library side univariate polynomials.
 * We have three types  of univariate polynomial depending on the plaintext space:
 * - PolyUnivZnx for coefficients in Z. We use signed 64-bits integers for the data.
 * - PolyUnivRnX for coefficients in R. We use double to approximate them;
 * - PolyUnivTnX for coefficients in T we scale from real in [0,1) to integers in [0, 2^64).
 * 
 * Only functions that prepares PolyUnivZnx are implemented, as both 
 * RnX and TnX univariate polynomial are converted to bivariate ZnX polynomial
 * before any operation.
 * 
 * For SVP-prepared univariate polynomials, it is used for the secret key which is
 * either a binary or a ternary key (-1,0,1), so we consider it to be a ZnX polynomial.
 */

// =============================================
// |                                           |
// |                  Types                    |
// |                                           |
// =============================================

/**
 * @brief Represents a \TnX polynomial, by means of
 * coefficients which are the numerators of a fraction with denominator \f$2^{64}\f$
 * @note The implicit length is \N for a polynomial in \TnX
 */
typedef uint64_t PolyUnivTnX;

/**
 * @brief Represents a \RnX polynomial.
 * @note The implicit length is \N for a polynomial in \RnX
 */
typedef double PolyUnivRnX;

/**
 * @brief Represents a \ZnX polynomial
 */
typedef int64_t PolyUnivZnx;

// =============================================
// |                                           |
// |                PolyUnivZnx                |
// |                                           |
// =============================================

/**
 * @brief Allocate a new univariate polynomial (in coefficient space).
 *
 * @param params_glwe The parameters to use for the polynomial.
 *
 * @retval NULL in case of a failure.
 * @retval A pointer to the new object otherwise.
 */
PolyUnivZnx* new_univ_znx(const GLWEParams* params_glwe);

/**
 * @brief Deallocates an univariate polynomial
 *
 * @param pol The object to deallocate
 */
void delete_univ_znx(PolyUnivZnx* pol);

// =============================================
// |                                           |
// |                    DFT                    |
// |                                           |
// =============================================

/**
 * @brief Allocate a new \ZnX univariate polynomial in DFT domain.
 *
 * @param module The backend instance.
 *
 * @retval NULL in case of a failure.
 * @retval An opaque type pointer to the new object otherwise.
 */
BE_PolyUnivZnxDFT* new_univ_znx_dft(const PvdaBackend* module);

/**
 * @brief Deallocates an \ZnX univariate polynomial in the DFT domain.
 *
 * @param pol The object to deallocate.
 */
void delete_univ_znx_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* pol);

/**
 * @brief Convert an univariate polynomial from coefficient space to DFT.
 *
 * @param module  The backend module.
 * @param res_dft The resulting DFT univariate polynomial.
 * @param a       The input univariate polynomial (in coefficient space).
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int univ_znx_to_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* res_dft, const PolyUnivZnx* a);

/**
 * @brief Convert an univariate polynomial from DFT domain into coefficient space.
 *
 * @param module The backend module.
 * @param res    The resulting coefficient-space univariate polynomial.
 * @param a_dft  The input DFT-domain polynomial.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int univ_dft_to_znx(const PvdaBackend* module, PolyUnivZnx* res, const BE_PolyUnivZnxDFT* a_dft);

// =============================================
// |                                           |
// |                   Prep                    |
// |                                           |
// =============================================

/**
 * @brief Allocate a new SVP-prepared univariate polynomial.
 *
 * @param module The backend instance.
 *
 * @retval NULL in case of a failure.
 * @retval An opaque type pointer to the new object otherwise.
 */
BE_PolyUnivSvpPrep* new_univ_svp_prep(const PvdaBackend* module);

/**
 * @brief Delete a SVP-prepared univariate polynomial.
 *
 * @param pol The object to deallocate.
 */
void delete_univ_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* pol);

/**
 * @brief Prepare an univariate polynomial.
 *
 * @param module   The backend module.
 * @param res_prep The resulting prepared (SVP) univariate polynomial.
 * @param a        The input univariate polynomial (in coefficient space).
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int univ_znx_to_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* res_prep, const PolyUnivZnx *a);

//* NOTE : We could have had a univ_prep_to_znx() function, but it is not mandatory.

// =============================================
// |                                           |
// |                   Tnx                     |
// |                                           |
// =============================================

/**
 * @brief Allocate a new Tnx univariate polynomial (in coefficient space).
 *
 * @param params_glwe The parameters to use for the polynomial.
 *
 * @retval NULL in case of a failure.
 * @retval A pointer to the new object otherwise.
 */
PolyUnivTnX* new_univ_tnx(const GLWEParams* params_glwe);

/**
 * @brief Deallocates a Tnx univariate polynomial.
 *
 * @param pol The object to deallocate.
 */
void delete_univ_tnx(PolyUnivTnX* pol);

// =============================================
// |                                           |
// |                   Rnx                     |
// |                                           |
// =============================================

/**
 * @brief Allocate a new Rnx univariate polynomial (in coefficient space).
 *
 * @param params_glwe The parameters to use for the polynomial.
 *
 * @retval NULL in case of a failure.
 * @retval a pointer to the new object otherwise.
 */
PolyUnivRnX* new_univ_rnx(const GLWEParams* params_glwe);

/**
 * @brief Deallocates a Rnx real univariate polynomial.
 *
 * @param pol The object to deallocate
 */
void delete_univ_rnx(PolyUnivRnX* pol);

// =============================================
// |                                           |
// |            Tnx-Rnx Conversions            |
// |                                           |
// =============================================

/**
 * @brief Converts a RnX polynomial into a TnX fixed-point one.
 *
 * @param params_glwe The GLWE parameters.
 * @param res         The resulting \TnX polynomial.
 * @param a           The input \RnX polynomial.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int univ_rnx_to_tnx(const GLWEParams* params_glwe, PolyUnivTnX* res, PolyUnivRnX* a);

/**
 * @brief Converts a TnX fixed-point polynomial into a RnX floating-point one.
 *
 * @param params_glwe The GLWE parameters.
 * @param res         The resulting \RnX polynomial.
 * @param a           The input \TnX polynomial.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int univ_tnx_to_rnx(const GLWEParams* params_glwe, PolyUnivRnX* res, PolyUnivTnX* a);

// =============================================
// |                                           |
// |                   Utils                   |
// |                                           |
// =============================================
/**
 * @brief Returns the number of bytes needed to store a univariate \ZnX polynomial.
 *
 * @param params_glwe The GLWE parameters.
 *
 * @return Number of bytes needed to store a univariate \ZnX polynomial with the given parameters.
 */
uint64_t univ_znx_bytes(const GLWEParams* params_glwe);

/**
 * @brief Returns the number of bytes needed to store a univariate real polynomial.
 *
 * @param params_glwe The GLWE parameters.
 *
 * @return Number of bytes needed to store a univariate \RnX polynomial with the given parameters.
 */
uint64_t univ_rnx_bytes(const GLWEParams* params_glwe);

/**
 * @brief Returns the number of bytes needed to store a univariate Tnx polynomial
 *
 * @param params_glwe The GLWE parameters.
 *
 * @return Number of bytes needed to store a univariate \TnX polynomial with the given parameters.
 */
uint64_t univ_tnx_bytes(const GLWEParams* params_glwe);

// =============================================
// |                                           |
// |              Print Functions              |
// |                                           |
// =============================================

/**
 * @brief Prints a Univariate Polynomial in \ZnX.
 *
 * @param pol A Pointer to the Univariate Polynomial.
 * @param nn  The size of the polynomial.
 */
void print_poly_univ_ZnX(PolyUnivZnx* pol, int64_t nn);

/**
 * @brief Prints a Univariate Polynomial in \RnX.
 *
 * @param pol A Pointer to the Univariate Polynomial.
 * @param nn  The size of the polynomial.
 */
void print_poly_univ_RnX(PolyUnivRnX* pol, int64_t nn);

/**
 * @brief Prints a Univariate Polynomial in \TnX.
 *
 * @param pol A Pointer to the Univariate Polynomial.
 * @param nn  The size of the polynomial.
 */
void print_poly_univ_TnX(PolyUnivTnX* pol, int64_t nn);

#endif
