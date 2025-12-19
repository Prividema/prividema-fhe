#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

/**
 * @file ggsw_ciphertext.h
 * @brief GGSW ciphertext
 * 
 * In this header file, we define the structure representing bivariate GGSW ciphertext 
 * and bivariate GGSW in DFT space. Moreover, we define alias for double and int64_t. 
 * Thus, we properly manipulate mathematical structure (polynomial, vector or matrix).
 * The types below are provided solely for clarity, nothing is changed in spqlios.
 */

#include <stdint.h>
#include <stdlib.h>
#include "ggsw_ct_params.h"

//! DEFINE SPQLIOS ALIAS (begin) 

/** @brief type that are included in */
typedef struct module_info_t MODULE;

/** @brief type that represents a prepared matrix */
typedef struct vmp_pmat_t VMP_PMAT;

/** @brief type that represents a bivariate polynomial in DFT space 
 *  @note `Length = l*N`
*/
typedef struct vec_znx_dft_t PolyBivDFT_;

/** @brief type that represents a bivariate polynomial vector in DFT space 
 *  @note `Length = (k+1) * l * N` for a bivGLWE
*/
typedef struct vec_znx_dft_t VecBivDFT_;

/** @brief type that represents a bivariate polynomial matrix in DFT space
 *  @note `Length = n_limbs_tilde * (k+1) * l * N` for a bivGGSW
*/
typedef struct vec_znx_dft_t MatBivDFT_;

/** @brief type that represents a univariate polynomial in DFT space 
 *  @note `Length = N for a GLWE`
*/
typedef struct vec_znx_dft_t PolyUnivDFT_;

/** @brief type that represents a univariate polynomial in DFT space 
 *  @note `Length = k * N for a GLWE`
*/
typedef struct vec_znx_dft_t VecUnivDFT_;


//! DEFINE GGSW ALIAS, NOT DFT & DFT (begin)

/** @brief type that represents a univariate polynomial.
 *  @note `Length = N for a polynomial in ZnX`.
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

/** @brief type that represents */

/** @brief type that represents a univariate polynomial in DFT space.
 *  @note `Length = N for a polynomial in ZnX`.
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


//! DEFINE bivGGSW AlIAS, NOT DFT & DFT (begin)

/** @brief type that represents a bivariate polynomial 
 *  @note `Length = l * N for bivariate polynomial in ZnXY`
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


typedef struct ggsw_ciphertext {
  GGSWCtParams* params;
  MatBiv* ct; //represent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in ZnX
} GGSWCiphertext;

int new_ggsw(GGSWCiphertext* res, GGSWCtParams* params, MatBiv* ct);
void delete_ggsw(GGSWCiphertext* ct);
void add_ggsw(GGSWCiphertext* res, GGSWCiphertext* ct1, GGSWCiphertext* ct2);
void const_mult_ggsw();


typedef struct ggsw_ciphertext_dft {
  GGSWCtParams* params;
  MatBivDFT* ct;
} GGSWPreparedCt;

int new_ggsw_prepared(GGSWPreparedCt* res, GGSWCtParams* params, MatBivDFT* ct);
void delete_ggsw_prepared(GGSWPreparedCt* ct);
void add_ggsw_prepared(GGSWCiphertext* res, GGSWCiphertext* ct1, GGSWCiphertext* ct2);

#endif  // GGSW_CIPHERTEXT_H
