#ifndef BIVARIATE_POLYNOMIAL_H
#define BIVARIATE_POLYNOMIAL_H

#include "glwe_params.h"
#include "maths_structures.h"
#include "spqlios_alias.h"

/**
 *
 * @file bivariate_polynomial.h
 *
 * @brief Operations on bivariate polynomials
 *
 */

// BIV POLY PART (begin)

/**
 * @brief Returns the number of coefficients in bivariate polynomial.
 *
 * @param params_glwe The bivariate GLWE parameters.
 * @return uint64_t
 *
 * @note The number of coefficients is the same in and out of the DFT domain.
 */
uint64_t poly_biv_coef_number(const GLWEParams* params_glwe);

/**
 * @brief Creates a bivariate polynomial
 *
 * @param params_glwe The bivGLWE parameters.
 */
PolyBiv* new_biv_poly(const GLWEParams* params_glwe);

void delete_biv(PolyBiv* biv);

/**
 * @brief Creates an allocated bivariate polynomial
 *
 * @param params_glwe The bivGLWE parameters.
 * @param biv_l The overriding l that will be used
 */
PolyBiv* new_biv_poly_custom_l(const GLWEParams* params_glwe, uint64_t biv_l);

/**
 * @brief Computes a random normal bivariate polynomial.
 *
 * By "random normal" it is meant that the coefficients of the polynomial before
 * the base-2k decomposition have been sampled from a normal distribution according
 * to the parameters
 *
 * @param result The result bivariate polynomial.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int normal_random_biv_poly(const GLWEParams* params_glwe, PolyBiv* result);

/**
 * @brief Computes a random uniform bivariate polynomial.
 *
 * By "random uniform" it is meant that the coefficients of the polynomial before
 * the base-2k decomposition have been sampled from a uniform distribution according
 * to the parameters
 *
 * @param params_glwe The bivGLWE parameters.
 * @param result     The output bivariate polynomial
 * @param precision The maximum degree in Y of the polynomial.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int uniform_random_biv_poly(const GLWEParams* params_glwe, PolyBiv* result, int64_t precision);

/**
 * @brief Adds two bivariate polynomials and puts the result in res.
 *
 * @param module The underlying compute module.
 * @param params_glwe The GLWE parameters.
 * @param res The result bivariate polynomial.
 * @param a The left-hand side bivariate polynomial.
 * @param b The right-hand side bivariate polynomial.
 */
void add_biv_poly(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBiv* a,
                  const PolyBiv* b);

/**
 * @brief Adds noise to the bivariate polynomial
 *
 * @param module The underlying compute module.
 * @param params_glwe The GLWE parameters (which decide the noise)
 * @param res The result bivariate polynomial, can be the same as the input
 * @param a The input bivariate polynomial.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int add_biv_noise(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBiv* a);

// BIV POLY IN DFT PART (begin)

/**
 * @brief The number of coefficients in bivariate polynomial.
 *
 * @param params_glwe The bivariate GLWE parameters.
 * @return The nubmer of coefficients in a bivariate polynomial
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * \ZnX, due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t poly_biv_coef_number_dft(const GLWEParams* params_glwe);

/**
 * @brief Creates an allocated bivariate polynomial in the DFT domain.
 *
 * @param params_glwe The bivGLWE parameters.
 */
PolyBivDFT* new_biv_poly_dft(const GLWEParams* params_glwe);

PolyBivDFT* new_biv_poly_dft_custom_l(const GLWEParams* params_glwe, uint64_t biv_l);

// COMMON PART (begin)

/**
 * @brief Returns the number of bytes needed to store a bivariate polynomial.
 *
 * @param params_glwe The GLWE parameters.
 * @return The bytes needed to store a bivariate polynomial with the given parameters
 *
 * @note The number of bytes needed to store a bivariate polynomial is the same in and out of the DFT domain.
 */
uint64_t poly_biv_bytes(const GLWEParams* params_glwe);

/**
 * @brief Computes P(X,2^(-kappa)) for P a bivariate polynomial.
 *
 * @param params_glwe The bivGLWE parameters.
 * @param res_univ The result univariate polynomial in \RnX.
 * @param pol The input bivariate polynomial.
 */
void biv_to_univ_rnx(const GLWEParams* params_glwe, PolyUnivRnX* res_univ, const PolyBiv* pol);

/**
 * @brief Computes the bivariate decompositionof a polynomial in \RnX.
 *
 * The output is quasi-normalized: the coefficients lie all in [-2^(K-1), 2^(K-1)] instead of
 * the normalized [-2^(K-1), 2^(K-1))
 *
 * @param params_glwe The bivGLWE parameters.
 * @param res The bivariate decomposition.
 * @param pol_univ The univariate polynomial.
 * @param bit_offset How many times the input should be right-shifted (divided by 2)
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int univ_rnx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivRnX* pol_univ, int64_t bit_offset);

/**
 * @brief Performs the DFT of a bivariate polynomial
 *
 *
 * @param module      The backend module
 * @param params_glwe The GLWE parameters for the associated bivariate polynomials
 * @param res_dft     The resulting DFT domain bivariate polynomial
 * @param a           The input coefficient-space bivariate polynomial
 *
 * @retval -1 If an error occurs
 * @retval 0  Otwerwise
 *
 */
int biv_coefs_to_dft(const MODULE* module, const GLWEParams* params_glwe, PolyBivDFT* res_dft, const PolyBiv* a);

/**
 * @brief Performs the iDFT of a bivariate polynomial
 *
 *
 * @param module      The backend module
 * @param params_glwe The GLWE parameters for the associated bivariate polynomials
 * @param res         The resulting coefficient-space bivariate polynomial
 * @param a_dft       The input DFT domain bivariate polynomial
 *
 * @retval -1 If an error occurs
 * @retval 0  Otwerwise
 *
 */
int biv_dft_to_coefs(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBivDFT* a_dft);

/**
 * @brief Computes P(X,2^(-kappa)) for P a bivariate polynomial. The result is in fixed-point representation.
 *
 * @param params_glwe The GLWE parameters.
 * @param res_tnx The result univariate polynomial in Tn[X].
 * @param pol The input bivariate polynomial.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int biv_to_univ_tnx(const GLWEParams* params_glwe, PolyUnivTnX* res_tnx, const PolyBiv* pol);

/**
 * @brief Computes the bivariate decomposition in Zn[X,Y] of a polynomial in Tn[X].
 *
 * The output is quasi-normalized: the coefficients lie all in [-2^(K-1), 2^(K-1)] instead of
 * the normalized [-2^(K-1), 2^(K-1))
 *
 * @param params_glwe The GLWE parameters.
 * @param res The bivariate decomposition.
 * @param pol_tnx The univariate polynomial in fixed-point form
 * @param bit_offset How many times the input should be right-shifted (divided by 2)
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int univ_tnx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivTnX* pol_tnx, int64_t bit_offset);

/**
 * @brief Computes the bivariate decomposition in Zn[X,Y] of a polynomial in Zn[X] divided by a power of 2.
 *
 * Since ZnX elements themselves don't fit in the Torus (or equivalently, they are all congruent to 0),
 * this function is only useful if called with bit_offset >= 1.
 *
 * The function is only well-defined for polynomials with coefficients in (-2^63, 2^63)
 * due to the internal conversion function used.
 *
 * The output is quasi-normalized: the coefficients lie all in [-2^(K-1), 2^(K-1)] instead of
 * the normalized [-2^(K-1), 2^(K-1))
 *
 * @param params_glwe The GLWE parameters.
 * @param res The bivariate decomposition.
 * @param pol_tnx The univariate ZnX polynomial with coefficients in (-2^62, 2^62)
 * @param bit_offset How many times the input should be right-shifted (divided by 2).
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int univ_znx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUniv* pol_univ, int64_t bit_offset);

#endif  // BIVARIATE_POLYNOMIAL_H
