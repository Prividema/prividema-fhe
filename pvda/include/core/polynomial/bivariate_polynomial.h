#ifndef BIVARIATE_POLYNOMIAL_H
#define BIVARIATE_POLYNOMIAL_H

#include "backend.h"
#include "glwe_params.h"
#include "pvda_types.h"
#include "univariate_polynomial.h"

/**
 * @file bivariate_polynomial.h
 *
 * @brief Operations on bivariate polynomials.
 * 
 * This file contains everything on library side bivariate polynomials.
 * 
 * We have a single type PolyBiv that is always a \Znx bivariate polynomial.
 * We convert Univariate polynomials over \ZnX, \RnX and \Tnx to PolyBiv with precision bounds.
 */

// =============================================
// |                                           |
// |                  PolyBiv                  |
// |                                           |
// =============================================
// The PolyBiv structure is defined in pvda_types.h

/**
 * @brief Creates a bivariate polynomial.
 *
 * @param params_glwe The GLWE parameters.
 */
PolyBiv* new_biv(const GLWEParams* params_glwe);

/**
 * @brief Creates a bivariate polynomial.
 *
 * @param nn Parameter N (number of coefficients per level).
 * @param biv_l Parameter l (number of levels / depth).
 */
PolyBiv* new_biv_custom_params(uint64_t nn, uint64_t biv_l);

/**
 * @brief Deletes a bivariate polynomial.
 *
 * @param biv The bivariate polynomial to delete.
 */
void delete_biv(PolyBiv* biv);

// =============================================
// |                                           |
// |                   Views                   |
// |                                           |
// =============================================
/**
 * @brief Creates a new bivariate polynomial view, that is, 
 *        a PolyBiv that does NOT own the memory it refers to.
 *
 * @param nn     The ring dimension.
 * @param l      The depth (parameter l) of the polynomial.
 * @param stride The distance in elements between limbs. For a contiguous layout,
 *               this is nn. In a GLWE, for example, it will usually be (k+1)*nn.
 * @param ptr    Pointer to the data of the polynomial.
 *
 * @return A PolyBiv view of the data with the given parameters. The return object
 *         should NOT have delete_biv called on.
 *
 */
PolyBiv biv_view(uint64_t nn, uint64_t l, int64_t stride, int64_t* ptr);


// =============================================
// |                                           |
// |                Operations                 |
// |                                           |
// =============================================

/** @brief res = a + b, limb by limb (no carry propagation). */
int biv_add(const PvdaBackend* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b);

/** @brief res = a - b, limb by limb (no carry propagation) */
int biv_sub(const PvdaBackend* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b);

/** @brief res = a(X^p), limb by limb */
int biv_automorphism(const PvdaBackend* module, int64_t p, PolyBiv* res, const PolyBiv* a);

/** @brief res = a with carry propagated across limbs in base 2^(log2_base)*/
int biv_normalize(const PvdaBackend* module, uint64_t log2_base, PolyBiv* res, const PolyBiv* a);

/**
 * @brief Adds noise to the bivariate polynomial.
 *
 * @param module      The underlying compute module.
 * @param params_glwe The GLWE parameters (which decide the noise).
 * @param res         The result bivariate polynomial, can be the same as the input.
 * @param a           The input bivariate polynomial.
 *
 * @note The noise is added to the last limb (limb l_a) only.
 * @note At present only NOISE_UNIFORM_POWER_OF_TWO is supported; other noise types return -1.
 * @note If res != a, both must be contiguous (stride == nn): out-of-place addition on strided
 *       polynomials is not implemented and returns -1.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int biv_add_noise(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBiv* a);

// =============================================
// |                                           |
// |                  Random                   |
// |                                           |
// =============================================

/**
 * @brief Computes a random normal bivariate polynomial.
 *
 * By "random normal" it is meant that the coefficients of the polynomial before
 * the base-2k decomposition have been sampled from a normal distribution according
 * to the parameters
 *
 * @param module The underlying compute module.
 * @param params_glwe The GLWE parameters.
 * @param res The result bivariate polynomial.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int biv_random_normal(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res);

/**
 * @brief Computes a random uniform bivariate polynomial.
 *
 * By "random uniform" it is meant that the coefficients of the polynomial before
 * the base-2k decomposition have been sampled from a uniform distribution according
 * to the parameters.
 *
 * @param module      The underlying compute module.
 * @param params_glwe The bivGLWE parameters.
 * @param res         The output bivariate polynomial.
 * @param precision   The number of limbs to fill (limbs 0 to precision-1), each sampled uniformly on K bits.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int biv_random_uniform(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res, int64_t precision);

// =============================================
// |                                           |
// |                Univ -> Biv                |
// |                                           |
// =============================================

/**
 * @brief Computes the bivariate decomposition of a polynomial in RnX.
 *
 * The output is quasi-normalized: the coefficients lie all in [-2^(K-1), 2^(K-1)] instead of
 * the normalized [-2^(K-1), 2^(K-1)).
 *
 * @param params_glwe The bivGLWE parameters.
 * @param res         The bivariate decomposition.
 * @param pol_univ    The univariate polynomial.
 * @param bit_offset  How many times the input should be right-shifted (divided by 2).
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int univ_rnx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivRnX* pol_univ, int64_t bit_offset);

/**
 * @brief Computes the bivariate decomposition in Zn[X,Y] of a polynomial in Tn[X].
 *
 * The output is quasi-normalized: the coefficients lie all in [-2^(K-1), 2^(K-1)] instead of
 * the normalized [-2^(K-1), 2^(K-1)).
 *
 * Note that the transformation internally works with 63 (instead of 64) bits of precision
 * for performance reasons.
 * The LSB of the provided TnX value is dropped (rounded down, which is not round-towards 0!).
 *
 *
 * @param params_glwe The GLWE parameters.
 * @param res         The bivariate decomposition.
 * @param pol_tnx     The univariate polynomial in fixed-point form.
 * @param bit_offset  How many times the input should be right-shifted (divided by 2).
 *                    Intended to be used both for GLWEGadget encryption as well as scaling factors for schemes.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int univ_tnx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivTnX* pol_tnx, int64_t bit_offset);

/**
 * @brief Computes the bivariate decomposition in Zn[X,Y] of a polynomial in Zn[X] divided by a power of 2.
 *
 * Since ZnX elements themselves don't fit in the Torus (or equivalently, they are all congruent to 0),
 * this function is only useful if called with bit_offset >= 1.
 *
 * The function is only well-defined for polynomials with coefficients in (-2^62, 2^62)
 * due to the internal conversion function used.
 *
 * The output is quasi-normalized: the coefficients lie all in [-2^(K-1), 2^(K-1)] instead of
 * the normalized [-2^(K-1), 2^(K-1)).
 *
 * @param params_glwe The GLWE parameters.
 * @param res         The bivariate decomposition.
 * @param pol_univ    The univariate ZnX polynomial with coefficients in (-2^62, 2^62).
 * @param bit_offset  How many times the input should be right-shifted (divided by 2).
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int univ_znx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivZnx* pol_univ, int64_t bit_offset);

/**
 * @brief Computes the bivariate decomposition in Zn[X,Y] of a multi-precision polynomial in Zn[X] divided by
 *        a power of 2, i.e. of a / 2^bit_offset mod 1.
 *
 * Unlike univ_znx_to_biv, the coefficients of a can be written over more than 64 bits. With
 * bit_offset = 64 * a->nb_words, a is read as a torus polynomial with denominator 2^(64 * nb_words)
 * (a multi-precision PolyUnivTnX).
 *
 * The value is rounded to the nearest multiple of 2^(-K * l_a) (ties rounded up). The output is normalized:
 * the coefficients of limbs 1..l_a-1 lie in [-2^(K-1), 2^(K-1)), and limb 0 too since the integer part is dropped.
 *
 * @param params_glwe The GLWE parameters (gives K and l_a).
 * @param res         The bivariate decomposition, with nn == a->nn and at least l_a limbs.
 *                    Only the first l_a limbs are written.
 * @param a           The multi-precision univariate polynomial.
 * @param bit_offset  How many times the input should be right-shifted (divided by 2). Can be negative.
 *
 * @retval -1 if the dimensions of res and a do not match.
 * @retval 0 otherwise.
 */
int univ_znx_mp_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivZnxMP* a, int64_t bit_offset);

// =============================================
// |                                           |
// |                Biv -> Univ                |
// |                                           |
// =============================================

/**
 * @brief Computes P(X,2^(-kappa)) for P a bivariate polynomial. The result is in fixed-point representation.
 *
 * @param params_glwe The GLWE parameters.
 * @param res_tnx     The result univariate polynomial in Tn[X].
 * @param pol         The input bivariate polynomial.
 * @param bit_offset  How many times the input has been right-shifted (divided by 2) / How many bits to left-shift (multiply by 2) the output.
 *                    Intended for use when the univ_to_biv has been used for a scaling factor.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int biv_to_univ_tnx(const GLWEParams* params_glwe, PolyUnivTnX* res_tnx, const PolyBiv* pol, int64_t bit_offset);

/**
 * @brief Computes P(X,2^(-kappa)) for P a bivariate polynomial.
 *
 * @param params_glwe The bivGLWE parameters.
 * @param res_univ    The result univariate polynomial in \RnX.
 * @param pol         The input bivariate polynomial.
 * @param bit_offset  How many times the input has been right-shifted (divided by 2) / How many bits to left-shift (multiply by 2) the output.
 *                    Intended for use when the univ_to_biv has been used for a scaling factor.
 */
void biv_to_univ_rnx(const GLWEParams* params_glwe, PolyUnivRnX* res_univ, const PolyBiv* pol, int64_t bit_offset);

/**
 * @brief Computes round(P(X, 2^(-K)) * 2^bit_offset) mod 2^(64 * res->nb_words) for P a bivariate polynomial,
 *        into a multi-precision \ZnX polynomial.
 *
 * This is the inverse of univ_znx_mp_to_biv: with bit_offset = 64 * res->nb_words, res is the multi-precision
 * torus representation of P. The conversion is exact (no rounding) when K * l_a <= bit_offset.
 *
 * @param params_glwe The GLWE parameters (gives K and l_a).
 * @param res         The resulting multi-precision polynomial, with nn == pol->nn.
 * @param pol         The input bivariate polynomial, with at least l_a limbs (only the first l_a are read).
 * @param bit_offset  How many bits to left-shift (multiply by 2) the output. Can be negative.
 *
 * @retval -1 if the dimensions of res and pol do not match or an allocation fails.
 * @retval 0 otherwise.
 */
int biv_to_univ_znx_mp(const GLWEParams* params_glwe, PolyUnivZnxMP* res, const PolyBiv* pol, int64_t bit_offset);

// =============================================
// |                                           |
// |                  Biv DFT                  |
// |                                           |
// =============================================

/**
 * @brief Allocates a bivariate polynomial in DFT domain (one DFT entry per limb, l_a limbs) through the backend.
 * @note The content is uninitialized.
 */
BE_PolyBivDFT* new_biv_dft(const PvdaBackend* module, const GLWEParams* params_glwe);

/**
 * @brief Allocates a bivariate polynomial in DFT domain with biv_l limbs through the backend.
 * @note The content is uninitialized.
 */
BE_PolyBivDFT* new_biv_dft_custom_params(const PvdaBackend* module, uint64_t biv_l);

/**
 * @brief Deletes a bivariate polynomial in DFT domain allocated by new_biv_dft*.
 */
void delete_biv_dft(const PvdaBackend* module, BE_PolyBivDFT* biv_dft);

// =============================================
// |                                           |
// |               Prepared Biv                |
// |                                           |
// =============================================

BE_PolyBivVmpPrep* new_biv_prep(const PvdaBackend* module, const GLWEParams* params_glwe);
BE_PolyBivVmpPrep* new_biv_prep_custom_params(const PvdaBackend* module, uint64_t biv_l);
void delete_biv_prep(const PvdaBackend* module, BE_PolyBivVmpPrep* biv_prep);

// =============================================
// |                                           |
// |                Transforms                 |
// |                                           |
// =============================================

/**
 * @brief Performs the DFT of a bivariate polynomial.
 *
 * @param module      The backend module.
 * @param params_glwe The GLWE parameters for the associated bivariate polynomials.
 * @param res_dft     The resulting DFT domain bivariate polynomial.
 * @param a           The input coefficient-space bivariate polynomial.
 *
 * @retval -1 If an error occurs.
 * @retval 0  Otherwise.
 *
 */
int biv_to_dft(const PvdaBackend* module, const GLWEParams* params_glwe, BE_PolyBivDFT* res_dft, const PolyBiv* a);

/**
 * @brief Performs the iDFT of a bivariate polynomial.
 *
 *
 * @param module      The backend module.
 * @param params_glwe The GLWE parameters for the associated bivariate polynomials.
 * @param res         The resulting coefficient-space bivariate polynomial.
 * @param a_dft       The input DFT domain bivariate polynomial.
 *
 * @retval -1 If an error occurs.
 * @retval 0  Otherwise.
 *
 */
int biv_from_dft(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res, const BE_PolyBivDFT* a_dft);

/**
 * @brief Transforms a bivariate polynomial into prepared form for operations like VMP.
 *
 * @param module      The backend module.
 * @param params_glwe The GLWE parameters for the associated bivariate polynomials.
 * @param res_prep    The resulting prepared (!= DFT) polynomial.
 * @param a           The input coefficient-space bivariate polynomial.
 *
 * @retval -1 If an error occurs.
 * @retval 0  Otherwise.
 */
int biv_to_prep(const PvdaBackend* module, const GLWEParams* params_glwe, BE_PolyBivVmpPrep* res_prep,
                      const PolyBiv* a);

// =============================================
// |                                           |
// |                   Utils                   |
// |                                           |
// =============================================

/**
 * @brief Returns the number of coefficients in bivariate polynomial.
 *
 * @param params_glwe The bivariate GLWE parameters.
 * @return The number of coefficients in bivariate polynomial.
 *
 * @note The number of coefficients is the same in and out of the DFT domain.
 */
uint64_t biv_coef_number(const GLWEParams* params_glwe);

/**
 * @brief The number of coefficients in bivariate polynomial.
 *
 * @param params_glwe The bivariate GLWE parameters.
 * @return The number of coefficients in a bivariate polynomial.
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * \ZnX, due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t biv_dft_coef_number(const GLWEParams* params_glwe);

/**
 * @brief Returns the number of bytes needed to store a bivariate polynomial.
 *
 * @param params_glwe The GLWE parameters.
 * @return The bytes needed to store a bivariate polynomial with the given parameters
 *
 * @note The number of bytes needed to store a bivariate polynomial is the same in and out of the DFT domain.
 */
uint64_t biv_bytes(const GLWEParams* params_glwe);

// =============================================
// |                                           |
// |                  Prints                   |
// |                                           |
// =============================================

/**
 * @brief Prints a Bivariate Polynomial.
 *
 * @param pol     A Pointer to the Bivariate Polynomial.
 */
void printf_poly_biv(PolyBiv* pol);

/**
 * @brief Prints a Bivariate Polynomial, with limits on the coefficients and limbs that will be printed
 *
 * @param pol     A Pointer to the Bivariate Polynomial.
 * @param max_n   The number of coefficients to print per limb
 * @param max_l   The number of limbs to print
 */
void print_coefs_biv(const PolyBiv* biv, int max_n, int max_l);


#endif  // BIVARIATE_POLYNOMIAL_H
