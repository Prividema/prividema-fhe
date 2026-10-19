#ifndef UNIVARIATE_POLYNOMIAL_H
#define UNIVARIATE_POLYNOMIAL_H

#include "backend.h"
#include "backend_types.h"
#include "glwe_params.h"
#include "maths_structures.h"

/**
 *
 * @file univariate_polynomial.h
 *
 * @brief Operations on univariate polynomials.
 *
 */

// =============================================
// |                                           |
// |                 PolyUniv                  |
// |                                           |
// =============================================

/**
 * @brief Allocate an new univariate polynomial (in coefficient space).
 *
 * @param params_glwe The parameters to use for the polynomial.
 *
 * @retval NULL in case of a failure.
 * @retval A pointer to the new object otherwise.
 */
PolyUniv* new_univ(const GLWEParams* params_glwe);

/**
 * @brief Deallocates an univariate polynomial
 *
 * @param pol The object to deallocate
 */
void delete_univ(PolyUniv* pol);

// =============================================
// |                                           |
// |                    DFT                    |
// |                                           |
// =============================================

/**
 * @brief Allocate a new univariate polynomial in DFT domain.
 *
 * @param module The backend instance.
 *
 * @retval NULL in case of a failure.
 * @retval An opaque type pointer to the new object otherwise.
 */
BE_PolyVecDFT* new_univ_dft(const PvdaBackend* module);

/**
 * @brief Deallocates an univariate polynomial in the DFT domain.
 *
 * @param pol The object to deallocate.
 */
void delete_univ_dft(const PvdaBackend* module, BE_PolyVecDFT* pol);

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
int univ_coefs_to_dft(const PvdaBackend* module, BE_PolyVecDFT* res_dft, const PolyUniv* a);

/**
 * @brief Convert an univariate polynomial from DFT domain into coefficient space.
 *
 * @param module The backend module.
 * @param res    The resulitng coefficient-space univariate polynomial.
 * @param a_dft  The input DFT-domain polynomial.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int univ_dft_to_coefs(const PvdaBackend* module, PolyUniv* res, const BE_PolyVecDFT* a_dft);

// =============================================
// |                                           |
// |                   Prep                    |
// |                                           |
// =============================================

/**
 * @brief Allocate a new prepared univariate polynomial.
 *
 * @param module The backend instance.
 *
 * @retval NULL in case of a failure.
 * @retval An opaque type pointer to the new object otherwise.
 */
BE_PolySvpPrep* new_univ_prep(const PvdaBackend* module);

/**
 * @brief Delete a prepared univariate polynomial.
 *
 * @param pol The object to deallocate.
 */
void delete_univ_prep(const PvdaBackend* module, BE_PolySvpPrep* pol);

/**
 * @brief Prepare an univariate polynomial.
 *
 * @param module   The backend module.
 * @param res_prep The resulitng coefficient-space univariate polynomial.
 * @param a        The input univariate polynomial (in coefficient space).
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int univ_coefs_to_prep(const PvdaBackend* module, BE_PolySvpPrep* res_prep, const PolyUniv *a);

//* NOTE : We could have had a univ_prep_to_coeff() function, but it is not mandatory.

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
 * @brief Converts a TnX polynomial into a RnX fixed-point one.
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
 * @return Number of bytes needed to store a univariate \ZnX polynomial wiht the given parameters.
 */
uint64_t poly_univ_bytes(const GLWEParams* params_glwe);

/**
 * @brief Returns the number of bytes needed to store a univariate real polynomial.
 *
 * @param params_glwe The GLWE parameters.
 *
 * @return Number of bytes needed to store a univariate \RnX polynomial wiht the given parameters.
 */
uint64_t poly_univ_rnx_bytes(const GLWEParams* params_glwe);

/**
 * @brief Returns the number of bytes needed to store a univariate Tnx polynomial
 *
 * @param params_glwe The GLWE parameters.
 *
 * @return Number of bytes needed to store a univariate \TnX polynomial wiht the given parameters.
 */
uint64_t poly_univ_tnx_bytes(const GLWEParams* params_glwe);

#endif
