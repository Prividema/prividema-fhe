#ifndef UNIVARIATE_POLYNOMIAL_H
#define UNIVARIATE_POLYNOMIAL_H

#include "glwe_params.h"
#include "maths_structures.h"
#include "spqlios_alias.h"

/*
 * @brief convert a univariate from coefficient space to DFT
 *
 * @param module The backend module
 * @param res_dft The resulting DFT domain univariate polynomial
 * @param a       The input univariate polynomial (in coefficient space)
 *
 */
int univ_coefs_to_dft(const MODULE* module, PolyUnivDFT* res_dft, const PolyUniv* a);

/**
 * @brief Convert a univariate polynomial from DFT domain into coefficient space
 *
 * @param module The backend module
 * @param res    The resulitng coefficient-space univariate polynomial
 * @param a_dft  The input DFT-domain polynomial
 */
int univ_dft_to_coefs(const MODULE* module, PolyUniv* res, const PolyUnivDFT* a_dft);

/**
 * @brief Returns the number of bytes needed to store a univariate polynomial.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return uint64_t
 *
 */
uint64_t poly_univ_bytes(const GLWEParams* params_glwe);

/**
 * @brief Returns the number of bytes needed to store a univariate real polynomial.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return uint64_t
 *
 */
uint64_t poly_univ_rnx_bytes(const GLWEParams* params_glwe);

/**
 * @brief Allocate a new univariate polynomial (in coefficient space)
 *
 * @param params_glwe The parameters to use for the polynomial
 *
 */
PolyUniv* new_univ(const GLWEParams* params_glwe);

/**
 * @brief Allocate a new univariate polynomial (in coefficient space)
 *
 * @param params_glwe The parameters to use for the polynomial
 *
 */
PolyUnivTnX* new_univ_tnx(const GLWEParams* params_glwe);

/**
 * @brief Allocate a new univariate polynomial in DFT domain
 *
 * @param params_glwe The parameters to use for the polynomial
 *
 */
PolyUnivDFT* new_univ_dft(const MODULE* module);

/**
 * @brief Allocate a new real univariate polynomial
 *
 * @param params_glwe The parameters to use for the polynomial
 *
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
void delete_univ_dft(PolyUnivDFT* pol);

int univ_rnx_to_tnx(const GLWEParams* params_glwe, PolyUnivTnX* res, PolyUnivRnX* a);

int univ_tnx_to_rnx(const GLWEParams* params_glwe, PolyUnivRnX* res, PolyUnivTnX* a);

#endif
