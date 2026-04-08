#include <stdint.h>

// =============================================
// |                                           |
// |         Variable presentation             |
// |                                           |
// =============================================

/**
 * In this library, the name of a variable follow this principle :
 *
 * var                : the elements of the variable are in the bivariate space.
 * var_dft            : the elements of the variable are bivariates elements in the DFT domain.
 * var_univariate     : the variable is a Zn[X] polynomial.
 * var_univariate_RnX : the variable is a Rn[X] polynomial.
 * var_univariate_dft : the variable is a Zn[X] polynomial in the DFT domain.
 */

// =============================================
// |                                           |
// |         Univariate Structures             |
// |                                           |
// =============================================

/**
 * @brief Represents a \f$\mathbb{T}_n[X]\f$ polynomial, by means of
 * coefficients which are the numerators of a fraction with denominator 2^64
 * @note `Length = N` for a polynomial in \f$\mathbb{T}_n[X]\f$.
 */
typedef uint64_t PolyUnivTnX;

/**
 * @brief Represents a vector of the above
 */
typedef uint64_t VecUnivTnX;

/**
 * @brief Represents a \f$\mathbb{R}_n[X]\f$ polynomial.
 * @note `Length = N` for a polynomial in \f$\mathbb{R}_n[X]\f$.
 */
typedef double PolyUnivRnX;

/**
 * @brief Represents a \f$\mathbb{R}_n[X]\f$ polynomial vector.
 * @note `Length = (k+1) * N` for a TRLWE.
 * @note `Length = k * N` for a TRLWE secret key.
 */
typedef double VecUnivRnX;

/** @brief Represents a \f$\mathbb{R}_n[X]\f$ polynomial vector.
 *  @note `Length = (k+1) * N` for a bivGLWE.
 *  @note `Length = k * N` for a bivGLWE.
 */
typedef int64_t VecUniv;

/** @brief Represents a \f$\mathbb{Z}_n[X]\f$ polynomial.
 *  @note `Length = N` for a polynomial in \f$\mathbb{Z}_n[X]\f$.
 */
typedef int64_t PolyUniv;

/** @brief Represents a \f$\mathbb{Z}_n[X]\f$ univariate polynomial vector.
 *  @note `Length = k * N` for a bivGLWE secret key.
 */
typedef int64_t VecUniv;

/** @brief Represents a univariate bivGGSW.
 *  @note `Length = n_limbs_tilde * k * N` for a bivGGSW.
 */
typedef int64_t MatUniv;

/** @brief Represents a \f$\mathbb{Z}_n[X]\f$ polynomial in the DFT domain.
 *  @note `Length = N` for a polynomial in \f$\mathbb{Z}_n[X]\f$.
 */
typedef double PolyUnivDFT;

/** @brief Represents a univariate \f$\mathbb{Z}_n[X]\f$ polynomial vector in the DFT domain.
 *  @note `Length = (k+1) * N` for a bivGLWE.
 *  @note `Length = k * N` for a bivGLWE secret key.
 */
typedef double VecUnivDFT;

/** @brief Represents a univariate bivGGSW in the DFT domain.
 *  @note `Length = n_limbs_tilde * (k+1) * N` for a bivGGSW.
 */
typedef double MatUnivDFT;

// =============================================
// |                                           |
// |          Bivariate Structures             |
// |                                           |
// =============================================

/** @brief Represents a bivariate polynomial.
 *  @note `Length = l * N` for bivariate polynomial in \f$\mathbb{Z}_n[X, Y]\f$.
 */
typedef int64_t PolyBiv;

/** @brief Represents a bivariate polynomial vector.
 *  @note `Length = l * (k+1) * N` for a bivGLWE.
 */
typedef int64_t VecBiv;

/** @brief Represents a bivariate polynomial matrix.
 *  @note `Length = n_limbs_tilde * l * (k+1) * N.
 */
typedef int64_t MatBiv;

/** @brief Represents a bivariate polynomial in the DFT domain.
 *  @note `Length = l*N`.
 */
typedef double PolyBivDFT;

/** @brief Represents a bivariate polynomial vector in the DFT domain.
 *  @note `Length = l * (k+1) * N` for a bivGLWE.
 */
typedef double VecBivDFT;

/** @brief Represents a bivariate polynomial matrix in the DFT domain.
 *  @note `Length = n_limbs_tilde * l * (k+1) * N` for a bivGGSW.
 */
typedef double MatBivDFT;

// =============================================
// |                                           |
// |              Print Functions              |
// |                                           |
// =============================================

/**
 * @brief Prints a Bivariate Polynomial.
 *
 * @param pol     A Pointer to the Bivariate Polynomial.
 * @param pol_sl  The stride length : The step to jump from the beggining to the next polynomial.
 * @param nn       The degree of the chosen cyclotomic polynomial.
 * @param l       The degree in Y.
 */
void printf_poly_biv(PolyBiv* pol, int64_t pol_sl, int64_t nn, int64_t l);

/**
 * @brief Prints a Vector of Bivariate Polynomial such as a bivGLWE ciphertext.
 *
 * @param pols          A Pointer to the vector.
 * @param pols_size     The size of the vector.
 * @param nn             The degree of the chosen cyclotomic polynomial.
 * @param l             The degree in Y.
 */
void printf_vec_poly_biv(VecBiv* pols, int64_t pols_size, int64_t nn, int64_t l);

/**
 * @brief Prints a Univariate Polynomial in \f$\mathbb{Z}_n[X]\f$.
 *
 * @param pol A Pointer to the Univariate Polynomial.
 * @param nn   The degree of the chosen cyclotomic polynomial.
 */
void printf_poly_univ_ZnX(PolyUniv* pol, int64_t nn);

/**
 * @brief Prints a Univariate Polynomial in \f$\mathbb{R}_n[X]\f$.
 *
 * @param pol A Pointer to the Univariate Polynomial.
 * @param nn   The degree of the chosen cyclotomic polynomial.
 */
void printf_poly_univ_RnX(double* pol, int64_t nn);

/**
 * @brief Prints a Vector of Univariate Polynomials.
 *
 * @param pols        A Pointer to the vector.
 * @param pols_size   The size of the vector.
 * @param nn           The degree of the chosen cyclotomic polynomial.
 */
void printf_vec_poly_univ(VecBiv* pols, int64_t pols_size, int64_t nn);

// TODO : Used for debugging. To be removed.
void printf_secret_key(PolyUniv** sk_values, int64_t nn, int64_t k);
