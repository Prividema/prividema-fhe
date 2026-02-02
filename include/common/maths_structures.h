#include <stdint.h>
//! DEFINE STRUCTURE ALIAS, NOT DFT & DFT (begin)

/** @brief type that represents a univariate polynomial.
 *  @note `Length = N for a polynomial in Zn[X]`.
*/
typedef int64_t PolyUniv;

/** @brief type that represents a univariate polynomial vector.
 *  @note `Length = k * N for a GLWE`.
*/
typedef int64_t VecUniv;

/** @brief type that represents a univariate GGSW.
 *  @note `Length = n_limbs_tilde * k * N for a GGSW`.
*/
typedef int64_t MatUniv;

/** @brief type that represents a univariate polynomial in DFT space.
 *  @note `Length = N for a polynomial in Zn[X]`.
*/
typedef double PolyUnivDFT;

/** @brief type that represents a univariate polynomial vector in DFT space.
 *  @note `Length = (k+1) * N for a GLWE`
*/
typedef double VecUnivDFT;

/** @brief type that represents a univariate GGSW in DFT space.
 *  @note `Length = n_limbs_tilde * (k+1) * N for a GGSW`
*/
typedef double MatUnivDFT;

//! DEFINE STRUCTURE ALIAS, NOT DFT & DFT (begin)

/** @brief type that represents a bivariate polynomial 
 *  @note `Length = l * N for bivariate polynomial in Zn[X,Y]`
*/
typedef int64_t PolyBiv;

/** @brief type that represents a bivariate polynomial vector 
 *  @note `Length = l * (k+1) * N for a bivGLWE`
*/
typedef int64_t VecBiv;

/** @brief type that represents a bivariate polynomial matrix 
 *  @note `Length = n_limbs_tilde * l * (k+1) * N
*/
typedef int64_t MatBiv;

/** @brief type that represents a bivariate polynomial in DFT space.
 *  @note `Length = l*N`.
*/
typedef double PolyBivDFT;

/** @brief type that represents a bivariate polynomial vector in DFT space.
 *  @note `Length = l * (k+1) * N` for a bivGLWE.
*/
typedef double VecBivDFT;

/** @brief type that represents a bivariate polynomial matrix in DFT space.
 *  @note `Length = n_limbs_tilde * l * (k+1) * N` for a bivGGSW.
*/
typedef double MatBivDFT;

//! DEFINE FUNCTION TO PRINT THEM

void printf_poly_biv(PolyBiv* pol, int64_t pol_sl, int64_t N, int64_t l);

void printf_vec_poly_biv(VecBiv* pols, int64_t pols_size, int64_t N, int64_t l);

void printf_poly_univ_ZnX(PolyUniv* pol, int64_t N);

void printf_poly_univ_RnX(double* pol, int64_t N);

void printf_vec_poly_univ(VecBiv* pols, int64_t pols_size, int64_t N);

void printf_secret_key(PolyUniv** sk_values, int64_t N, int64_t k);