#ifndef UNIVARIATE_POLYNOMIAL_H
#define UNIVARIATE_POLYNOMIAL_H

#include "backend.h"
#include "glwe_params.h"
#include "maths_structures.h"

/**
 *
 * @ file univariate_polynomial.h
 *
 * @brief Oerations on univariate polynomials
 *
 */

/*
 * @brief convert a univariate polynomial from coefficient space to DFT
 *
 * @param module The backend module
 * @param res_dft The resulting DFT univariate polynomial
 * @param a       The input univariate polynomial (in coefficient space)
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int univ_coefs_to_dft(const PvdaBackend* module, PolyUnivDFT* res_dft, const PolyUniv* a);

/**
 * @brief Convert a univariate polynomial from DFT domain into coefficient space
 *
 * @param module The backend module
 * @param res    The resulitng coefficient-space univariate polynomial
 * @param a_dft  The input DFT-domain polynomial
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int univ_dft_to_coefs(const PvdaBackend* module, PolyUniv* res, const PolyUnivDFT* a_dft);

/**
 * @brief Returns the number of bytes needed to store a univariate \ZnX polynomial.
 *
 * @param params_glwe The GLWE parameters.
 *
 * @return Number of bytes needed to store a univariate \ZnX polynomial wiht the given parameters
 */
uint64_t poly_univ_bytes(const GLWEParams* params_glwe);

/**
 * @brief Returns the number of bytes needed to store a univariate real polynomial.
 *
 * @param params_glwe The GLWE parameters.
 *
 * @return Number of bytes needed to store a univariate \RnX polynomial wiht the given parameters
 */
uint64_t poly_univ_rnx_bytes(const GLWEParams* params_glwe);

/**
 * @brief Allocate a new univariate polynomial (in coefficient space)
 *
 * @param params_glwe The parameters to use for the polynomial
 *
 * @return NULL in case of a failure, a pointer to the new object otherwise
 */
PolyUniv* new_univ(const GLWEParams* params_glwe);

/**
 * @brief Allocate a new univariate polynomial (in coefficient space)
 *
 * @param params_glwe The parameters to use for the polynomial
 *
 * @return NULL in case of a failure, a pointer to the new object otherwise
 */
PolyUnivTnX* new_univ_tnx(const GLWEParams* params_glwe);

/**
 * @brief Allocate a new univariate polynomial in DFT domain
 *
 * @param module The underlying compute module where the new DFT univariate polynomial should be allocated
 *
 * @return NULL in case of a failure, a pointer to the new object otherwise
 */
PolyUnivDFT* new_univ_dft(const PvdaBackend* module);

/**
 * @brief Allocate a new real univariate polynomial
 *
 * @param params_glwe The parameters to use for the polynomial
 *
 * @return NULL in case of a failure, a pointer to the new object otherwise
 */
PolyUnivRnX* new_univ_rnx(const GLWEParams* params_glwe);

/**
 * @brief Deallocates a univariate polynomial
 *
 * @param pol The object to deallocate
 */
void delete_univ(PolyUniv* pol);

/**
 * @brief Deallocates a univariate polynomial
 *
 * @param pol The object to deallocate
 */
void delete_univ_tnx(PolyUnivTnX* pol);

/**
 * @brief Deallocates a real univariate polynomial
 *
 * @param pol The object to deallocate
 */
void delete_univ_rnx(PolyUnivRnX* pol);

/**
 * @brief Deallocates a univariate polynomial in the DFT domain
 *
 * @param pol The object to deallocate
 */
void delete_univ_dft(const PvdaBackend* module, PolyUnivDFT* pol);

/**
 * @brief Converts a \RnX polynomial into a \TnX fixed-point one
 *
 * @param params_glwe The GLWE parameters
 * @param res         The resulting \TnX polynomial
 * @param a           The input \RnX polynomial
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int univ_rnx_to_tnx(const GLWEParams* params_glwe, PolyUnivTnX* res, PolyUnivRnX* a);

/**
 * @brief Converts a \TnX polynomial into a \RnX fixed-point one
 *
 * @param params_glwe The GLWE parameters
 * @param res         The resulting \RnX polynomial
 * @param a           The input \TnX polynomial
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int univ_tnx_to_rnx(const GLWEParams* params_glwe, PolyUnivRnX* res, PolyUnivTnX* a);

#endif
