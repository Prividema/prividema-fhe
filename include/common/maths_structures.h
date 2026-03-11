#include <stdint.h>

// =============================================
// |                                           |
// |         Univariate Structures             |
// |                                           |
// =============================================

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
 *  @note `Length = (k+1) * N` for a GLWE.
 *  @note `Length = k * N` for a GLWE.
 */
typedef int64_t VecUniv;

/** @brief Represents a \f$\mathbb{Z}_n[X]\f$ polynomial.
 *  @note `Length = N` for a polynomial in \f$\mathbb{Z}_n[X]\f$.
 */
typedef int64_t PolyUniv;

/** @brief Represents a \f$\mathbb{Z}_n[X]\f$ univariate polynomial vector.
 *  @note `Length = k * N` for a GLWE secret key.
 */
typedef int64_t VecUniv;

/** @brief Represents a univariate GGSW.
 *  @note `Length = n_limbs_tilde * k * N` for a GGSW.
 */
typedef int64_t MatUniv;

/** @brief Represents a \f$\mathbb{Z}_n[X]\f$ polynomial in the DFT domain.
 *  @note `Length = N` for a polynomial in \f$\mathbb{Z}_n[X]\f$.
 */
typedef double PolyUnivDFT;

/** @brief Represents a univariate \f$\mathbb{Z}_n[X]\f$ polynomial vector in the DFT domain.
 *  @note `Length = (k+1) * N` for a GLWE.
 *  @note `Length = k * N` for a GLWE secret key.
 */
typedef double VecUnivDFT;

/** @brief Represents a univariate GGSW in the DFT domain.
 *  @note `Length = n_limbs_tilde * (k+1) * N` for a GGSW.
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
 * @param N       The degree of the chosen cyclotomic polynomial.
 * @param l       The degree in Y.
 */
void printf_poly_biv(PolyBiv* pol, int64_t pol_sl, int64_t N, int64_t l);

/**
 * @brief Prints a Vector of Bivariate Polynomial such as a GLWE ciphertext.
 * 
 * @param pols          A Pointer to the vector.
 * @param pols_size     The size of the vector.
 * @param N             The degree of the chosen cyclotomic polynomial.
 * @param l             The degree in Y.
 */
void printf_vec_poly_biv(VecBiv* pols, int64_t pols_size, int64_t N, int64_t l);

/**
 * @brief Prints a Univariate Polynomial in \f$\mathbb{Z}_n[X]\f$.
 * 
 * @param pol A Pointer to the Univariate Polynomial.
 * @param N   The degree of the chosen cyclotomic polynomial.
 */
void printf_poly_univ_ZnX(PolyUniv* pol, int64_t N);

/**
 * @brief Prints a Univariate Polynomial in \f$\mathbb{R}_n[X]\f$.
 * 
 * @param pol A Pointer to the Univariate Polynomial.
 * @param N   The degree of the chosen cyclotomic polynomial.
 */
void printf_poly_univ_RnX(double* pol, int64_t N);

/**
 * @brief Prints a Vector of Univariate Polynomials.
 * 
 * @param pols        A Pointer to the vector.
 * @param pols_size   The size of the vector.
 * @param N           The degree of the chosen cyclotomic polynomial.
 */
void printf_vec_poly_univ(VecBiv* pols, int64_t pols_size, int64_t N);

// TODO : Used for debugging. To be removed.
void printf_secret_key(PolyUniv** sk_values, int64_t N, int64_t k);