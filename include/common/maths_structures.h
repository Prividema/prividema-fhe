#include <stdint.h>

// =============================================
// |                                           |
// |         Variable presentation             |
// |                                           |
// =============================================

/*
 * In this library, an attempt is made to name variables as follows :
 *
 * var          : the elements of the variable are in the bivariate space.
 * var_dft      : the elements of the variable are bivariates elements in the DFT domain.
 * var_univ     : the variable is a \ZnX polynomial.
 * var_univ_rnx : the variable is a \RnX polynomial.
 * var_univ_tnx : the variable is a \TnX polynomial using fixed point representation.
 * var_univ_dft : the variable is a \ZnX polynomial in the DFT domain.
 */

// =============================================
// |                                           |
// |         Univariate Structures             |
// |                                           |
// =============================================

/**
 * @brief Represents a \TnX polynomial, by means of
 * coefficients which are the numerators of a fraction with denominator \f$2^{64}\f$
 * @note The implicit lenght is \N for a polynomial in \TnX
 */
typedef uint64_t PolyUnivTnX;

/**
 * @brief Represents a vector of univariate \TnX polynomials (flattened, each encoded as in PolyUnivTnX)
 */
typedef uint64_t VecUnivTnX;

/**
 * @brief Represents a \RnX polynomial.
 * @note The implicit lenght is \N for a polynomial in \RnX
 */
typedef double PolyUnivRnX;

/**
 * @brief Represents a vector of univariate \RnX polynomial (flattened, each encoded as in PolyUnivRnX)
 */
typedef double VecUnivRnX;

/**
 * @brief Represents a \ZnX polynomial
 */
typedef int64_t PolyUniv;

/**
 * @brief Represents a vector of \ZnX polynomials (flattened, each encoded as in PolyUniv)
 */
typedef int64_t VecUniv;

/**
 * @brief Represents a univariate polynomial in the DFT domain.
 */
typedef double PolyUnivDFT;

/**
 * @brief Represents a vector of univariate polynomials in the DFT domain (flattened, each encoded as in PolyUnivDFT)
 */
typedef double VecUnivDFT;

// =============================================
// |                                           |
// |          Bivariate Structures             |
// |                                           |
// =============================================

/**
 * @brief Represents a bivariate polynomial.
 */
typedef int64_t PolyBiv;

/**
 * @brief Represents a bivariate polynomial vector (flattened).
 */
typedef int64_t VecBiv;

/**
 * @brief Represents a bivariate polynomial matrix (flattened)
 */
typedef int64_t MatBiv;

/**
 * @brief Represents a bivariate polynomial in the DFT domain.
 */
typedef double PolyBivDFT;

/**
 *  @brief Represents a bivariate polynomial vector in the DFT domain.
 */
typedef double VecBivDFT;

/**
 * @brief Represents a bivariate polynomial matrix in the DFT domain.
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
 * @param nn      The degree of the chosen cyclotomic polynomial.
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
