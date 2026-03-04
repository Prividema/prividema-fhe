#ifndef BIVARIATE_POLYNOMIAL_H
#define BIVARIATE_POLYNOMIAL_H

#include <stdint.h>

#include "glwe_ct_params.h"
#include "maths_structures.h"
#include "spqlios_alias.h"

//! BIV POLY PART (begin)

/**
 * @brief Returns the number of coefficient in bivariate polynomial.
 *
 * @param params The GLWE parameters.
 * @return int64_t
 *
 * @note The number of coefficient is the same in and out DFT space.
 */
uint64_t poly_biv_coef_number(const GLWECtParams* params);

/**
 * @brief Computes a random normal bivariate polynomial.
 *
 * @param module The module holding the degree N and FFT64.
 * @param params The GLWE parameters.
 * @return PolyBiv*
 */
PolyBiv* new_normal_random_biv_poly(const MODULE* module, const GLWECtParams* params);

/**
 * @brief Computes a random uniform bivariate polynomial.
 *
 * @param module The module holding the degree N and FFT64.
 * @param params The GLWE parameters.
 * @param precision The maximum degree in Y of the polynomial.
 * @return PolyBiv*
 */
PolyBiv* new_uniform_random_biv_poly(const MODULE* module, const GLWECtParams* params, int64_t precision);

/**
 * @brief Adds two bivariate polynomial and puts it in res.
 *
 * @param params The GLWE parameters.
 * @param res The result bivariate polynomial.
 * @param res_sl The stride between each Zn[X] polynomial in res.
 * @param a The left-hand side bivariate polynomial.
 * @param a_sl The stride between each Zn[X] polynomial in a.
 * @param b The right-hand side bivariate polynomial.
 * @param b_sl The stride between each Zn[X] polynomial in b.
 */
void add_biv_poly(const GLWECtParams* params, PolyBiv* res, int64_t res_sl, const PolyBiv* a, int64_t a_sl, const PolyBiv* b,
                  int64_t b_sl);

//! BIV POLY IN DFT PART (begin)

/**
 * @brief The number of coefficient in bivariate polynomial.
 *
 * @param params The GLWE parameters.
 * @return int64_t
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * Zn[X], due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t poly_biv_coef_number_dft(const GLWECtParams* params);

/**
 * @brief Computes a random normal bivariate polynomial in the DFT domain.
 *
 * @param module The module holding the degree N and FFT64.
 * @param params The GLWE parameters.
 * @return PolyBivDFT*
 */
PolyBivDFT* new_normal_random_biv_poly_dft(const MODULE* module, const GLWECtParams* params);

/**
 * @brief Computes a random uniform bivariate polynomial.
 *
 * @param module The module holding the degree N and FFT64.
 * @param params The GLWE parameters.
 * @return PolyBiv*
 */
PolyBivDFT* new_uniform_random_biv_poly_dft(const MODULE* module, const GLWECtParams* params, int64_t precision);

/**
 * @brief Adds two bivariate polynomial and puts it in res in the DFT domain.
 *
 * @param params The GLWE parameters.
 * @param res_dft The result bivariate polynomial in the DFT domain.
 * @param res_sl The stride between each Zn[X] polynomialin res_dft.
 * @param a_dft The left-hand side bivariate polynomial in the DFT domain.
 * @param a_sl The stride between each Zn[X] polynomial in a_dft.
 * @param b_dft The right-hand side bivariate polynomial in the DFT domain.
 * @param b_sl The stride between each Zn[X] polynomial in b_dft.
 */
void add_biv_poly_dft(const GLWECtParams* params, PolyBivDFT* res_dft, int64_t res_sl, const PolyBivDFT* a_dft, int64_t a_sl,
                      const PolyBivDFT* b_dft, int64_t b_sl);

//! COMMON PART (begin)

/**
 * @brief Returns the number of bytes needed to store a bivariate polynomial.
 *
 * @param params The GLWE parameters.
 * @return uint64_t
 *
 * @note The number of bytes needed to store a bivariate polynomial is the same in and out of DFT space.
 */
uint64_t poly_biv_bytes(const GLWECtParams* params);

// TODO modifie size partout
/**
 *
 * @brief Returns the size - the number l of Zn[X] coefficients - for a bivariate polynomial.
 *
 * @param params
 * @return uint64_t
 *
 * @note The size of a bivariate polynomial is the same in and out of DFT space.
 */
uint64_t poly_biv_size(const GLWECtParams* params);

/**
 * @brief Returns the number of bytes needed to store a univariate polynomial.
 *
 * @param params The GLWE parameters.
 * @return uint64_t
 *
 * @note The number of bytes needed to store an univariate polynomial is the same in and out of DFT space.
 */
uint64_t poly_univ_bytes(const GLWECtParams* params);

/**
 * @brief Computes P(X,2^(-kappa)) for P a bivariate polynomial.
 *
 * @param params The GLWE parameters.
 * @param pol_univ The result univariate polynomial in Rn[X].
 * @param pol_biv The input bivariate polynomial.
 */
void biv_to_univ(const GLWECtParams* params, double* res_univ, const PolyBiv* pol);

/**
 * @brief Computes the bivariate decomposition in Zn[X,Y] of a polynomial in Rn[X].
 *
 * @param params The GLWE parameters.
 * @param res The bivariate decomposition.
 * @param pol_univ The univariate polynomial.
 * @return int
 */
int univ_to_biv(const GLWECtParams* params, PolyBiv* res, const double* pol_univ);

#endif  // BIVARIATE_POLYNOMIAL_H