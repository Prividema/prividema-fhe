#ifndef MATHS_STRUCTURES_H
#define MATHS_STRUCTURES_H

#include <stdint.h>

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
 * @brief Represents a \RnX polynomial.
 * @note The implicit lenght is \N for a polynomial in \RnX
 */
typedef double PolyUnivRnX;

/**
 * @brief Represents a \ZnX polynomial
 */
typedef int64_t PolyUniv;



// =============================================
// |                                           |
// |          Bivariate Structures             |
// |                                           |
// =============================================

/**
 * @brief Structure for bivariate polynomials.
 * 
 * Owners are created with new_biv/new_biv_custom_params and freed with delete_biv.
 * 
 * Views are created with biv_view ot the *_extract_poly_view functions and must not be deleted.
 */
typedef struct poly_biv_t
{
	uint64_t nn;       	///< Ring dimension (also number of coefficient in one polynomial)
	uint64_t nb_limbs; 	///< The number of limbs (base-2^K digits) of the polynomial
	int64_t stride;  	///< Stride (distance) between consecutive limbs of the bivariate polynomial, in elements.
	                 	///< For a contiguous layout, this will be nn. For a GLWE with k=2, for example, it would be 2*nn
	int64_t* ptr;  		///< Pointer to the underlying data
} PolyBiv;

// =============================================
// |                                           |
// |              Print Functions              |
// |                                           |
// =============================================

/**
 * @brief Prints a Bivariate Polynomial.
 *
 * @param pol     A Pointer to the Bivariate Polynomial.
 */
void printf_poly_biv(PolyBiv* pol);

/**
 * @brief Prints a Bivariate Polynomial, with limits on the coefifcients and limbs that will be printed
 *
 * @param pol     A Pointer to the Bivariate Polynomial.
 * @param max_n   The number of coefficients to print per limb
 * @param max_l   The number of limbs to print
 */
void print_coefs_biv(const PolyBiv* biv, int max_n, int max_l);

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
void printf_vec_poly_univ(int64_t* pols, int64_t pols_size, int64_t nn);

#endif
