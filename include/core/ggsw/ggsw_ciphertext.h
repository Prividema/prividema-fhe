#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

#include <stdint.h>
#include <stdlib.h>
#include "ggsw_ct_params.h"

// We are not defining types that differ from SPQLios. 
// The types below are provided solely for clarity. 
// If there is "_" at the end, we just choose a new alias for an opaque type.

/** DEFINE SPQLIOS ALIAS (begin) */

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

/** DEFINE SPQLIOS ALIAS (end) */ 


/** DEFINE GGSW ALIAS, NOT DFT & DFT (begin) */

/** @brief type that represents a univariate polynomial 
 *  @note `Length = N for a polynomial in ZnX`
*/
typedef int64_t PolyUniv;

/** @brief type that represents a bivariate polynomial vector 
 *  @note `Length = k * N for a GLWE`
*/
typedef int64_t VecUniv;

/** @brief type that represents */

/** @brief type that represents a univariate polynomial in DFT space 
 *  @note `Length = N for a polynomial in ZnX`
*/
typedef double PolyUnivDFT;

/** @brief type that represents a univariate polynomial in DFT space 
 *  @note `Length = k * N for a GLWE`
*/
typedef double VecUnivDFT;

/** DEFINE GGSW ALIAS, NOT DFT & DFT (end) */

/** @brief type that represents a bivariate polynomial in DFT space 
 *  @note `Length = l*N`
*/
typedef double PolyBivDFT;

/** @brief type that represents a bivariate polynomial vector in DFT space 
 *  @note `Length = (k+1) * l * N` for a bivGLWE
*/
typedef double VecBivDFT;

/** @brief type that represents a bivariate polynomial matrix in DFT space
 *  @note `Length = n_limbs_tilde * (k+1) * l * N` for a bivGGSW
*/
typedef double MatBivDFT;



/** @brief type that represents a bivariate polynomial 
 *  @note `Length = l * N for bivariate polynomial in ZnXY`
*/
typedef int64_t PolyBiv;

/** @brief type that represents a bivariate polynomial vector 
 *  @note `Length = l * k * N for a bivGLWE`
*/
typedef int64_t VecBiv;

/** @brief type that represents a bivariate polynomial matrix 
 *  @note `Length = n_limbs_tilde * l * k * N
*/
typedef int64_t MatBiv;



// GGSW encode small integer polynomial represented via a function f
 // We represent GGSW as a matrix of size n_limbs x n_limbs_tilde with coefficients that are in ZnX
 // each row i is a GLWE ciphertext that encrypts f(H(i)).
 // [ GLWE(f(H(0,l_tilde)))
 // ...
 // GLWE(f(H(n_limbs_tilde-1,l_tilde)))]

 // H is the gadget matrix (bg=2^{kappa_tilde}) reorganized in power of bg (in order to keep the prefix property)
 // [1/bg 0 0 ...0 0] 
 // [0 1/bg 0 ...0 0]
 // ...
 // [0 0 0... 0 1/bg]
 // [1/bg² 0 0...0 0]
 // [0 1/bg² 0 ..0 0]
 // ...
 // [0 0 ... 0 1/bg²]
 //...
 // [1/bg^l_t 0 ...0 0]
 //...
 // [0 0... 0 1/bg^l_t]

typedef struct ggsw_ciphertext {
  GGSWCtParams* params;
  MatBiv* ct; //represent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in ZnX
} GGSWCiphertext;

int new_ggsw(GGSWCiphertext* res, GGSWCtParams* params, MatBiv* ct);
void delete_ggsw(GGSWCiphertext* ct);

typedef struct ggsw_ciphertext_dft {
  GGSWCtParams* params;
  MatBivDFT* ct;
} GGSWPreparedCt;

int new_ggsw_prepared(GGSWPreparedCt* res, GGSWCtParams* params, MatBivDFT* ct);
void delete_ggsw_prepared(GGSWPreparedCt* ct);

#endif  // GGSW_CIPHERTEXT_H
