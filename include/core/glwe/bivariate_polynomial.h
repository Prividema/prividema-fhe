#ifndef BIVARIATE_POLYNOMIAL_H
#define BIVARIATE_POLYNOMIAL_H

#include "glwe_params.h"
#include "spqlios_alias.h"

//! BIV POLY PART (begin)

/**
 * @brief Returns the number of coefficient in bivariate polynomial.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return uint64_t
 *
 * @note The number of coefficient is the same in and out of the DFT domain.
 */
uint64_t poly_biv_coef_number(const GLWEParams* params_glwe);

/**
 * @brief Creates an allocated bivariate polynomial
 *
 * @param params_glwe The bivGLWE parameters.
 */
PolyBiv* new_biv_poly(const GLWEParams* params_glwe);

/**
 * @brief Computes a random normal bivariate polynomial.
 *
 * By "random normal" it is meant that the coefficients of the polynomial before
 * the base-2k decomposition have been sampled from a normal distribution according
 * to the parameters
 *
 * @param params_glwe The bivGLWE parameters.
 * @param result The result bivariate polynomial.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int normal_random_biv_poly(const GLWEParams* params_glwe, PolyBiv* result);

/**
 * @brief Computes a random uniform bivariate polynomial.
 *
 * By "random uniform" it is meant that the coefficients of the polynomial before
 * the base-2k decomposition have been sampled from a distribution according
 * to the parameters
 *
 * @param params_glwe The bivGLWE parameters.
 * @param precision The maximum degree in Y of the polynomial.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int uniform_random_biv_poly(const GLWEParams* params_glwe, PolyBiv* result, int64_t precision);

/**
 * @brief Adds two bivariate polynomial and puts it in res.
 *
 * @param params_glwe The bivGLWE parameters.
 * @param res The result bivariate polynomial.
 * @param res_sl The stride between each Zn[X] polynomial in res.
 * @param a The left-hand side bivariate polynomial.
 * @param a_sl The stride between each Zn[X] polynomial in a.
 * @param b The right-hand side bivariate polynomial.
 * @param b_sl The stride between each Zn[X] polynomial in b.
 */
void add_biv_poly(const GLWEParams* params_glwe, PolyBiv* res, int64_t res_sl, const PolyBiv* a, int64_t a_sl,
                  const PolyBiv* b, int64_t b_sl);

//! BIV POLY IN DFT PART (begin)

/**
 * @brief The number of coefficient in bivariate polynomial.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return uint64_t
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * Zn[X], due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t poly_biv_coef_number_dft(const GLWEParams* params_glwe);

/**
 * @brief Creates an allocated bivariate polynomial in the DFT domain.
 *
 * @param params_glwe The bivGLWE parameters.
 */
PolyBivDFT* new_biv_poly_dft(const GLWEParams* params_glwe);

/**
 * @brief Computes a random normal bivariate polynomial in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param params_glwe The bivGLWE parameters.
 * @param result_dft The result bivarariate polynomial in the DFT domain.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int normal_random_biv_poly_dft(const MODULE* module, const GLWEParams* params_glwe, PolyBivDFT* result_dft);

/**
 * @brief Computes a random uniform bivariate polynomial.
 *
 * @param module Additionnal information for backend.
 * @param params_glwe The bivGLWE parameters.
 * @param result_dft The result bivariate polynomial in the DFT domain.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int uniform_random_biv_poly_dft(const MODULE* module, const GLWEParams* params_glwe, PolyBivDFT* result_dft,
                                int64_t precision);

/**
 * @brief Adds two bivariate polynomial and puts it in res in the DFT domain.
 *
 * @param params_glwe The bivGLWE parameters.
 * @param res_dft The result bivariate polynomial in the DFT domain.
 * @param res_sl The stride between each Zn[X] polynomialin res_dft.
 * @param a_dft The left-hand side bivariate polynomial in the DFT domain.
 * @param a_sl The stride between each Zn[X] polynomial in a_dft.
 * @param b_dft The right-hand side bivariate polynomial in the DFT domain.
 * @param b_sl The stride between each Zn[X] polynomial in b_dft.
 */
void add_biv_poly_dft(const GLWEParams* params_glwe, PolyBivDFT* res_dft, int64_t res_sl, const PolyBivDFT* a_dft,
                      int64_t a_sl, const PolyBivDFT* b_dft, int64_t b_sl);

//! COMMON PART (begin)

/**
 * @brief Returns the number of bytes needed to store a bivariate polynomial.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return uint64_t
 *
 * @note The number of bytes needed to store a bivariate polynomial is the same in and out of the DFT domain.
 */
uint64_t poly_biv_bytes(const GLWEParams* params_glwe);

/**
 * @brief Returns the size - the number l of Zn[X] coefficients - for a bivariate polynomial.
 *
 * @param params_glwe
 * @return uint64_t
 *
 * @note The size of a bivariate polynomial is the same in and out of the DFT domain.
 */
uint64_t poly_biv_size(const GLWEParams* params_glwe);

/**
 * @brief Computes P(X,2^(-kappa)) for P a bivariate polynomial.
 *
 * @param params_glwe The bivGLWE parameters.
 * @param pol_univ The result univariate polynomial in Rn[X].
 * @param pol_biv The input bivariate polynomial.
 */
void biv_to_univ(const GLWEParams* params_glwe, double* res_univ, const PolyBiv* pol);

/**
 * @brief Computes the bivariate decomposition in Zn[X,Y] of a polynomial in Rn[X].
 *
 * @param params_glwe The bivGLWE parameters.
 * @param res The bivariate decomposition.
 * @param pol_univ The univariate polynomial.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int univ_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const double* pol_univ);

int biv_coefs_to_dft(const MODULE* module, const GLWEParams* params_glwe, PolyBivDFT* res_dft, const PolyBiv* a);

int biv_dft_to_coefs(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBivDFT* a_dft);

#endif  // BIVARIATE_POLYNOMIAL_H
