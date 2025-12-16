#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

#include <stdint.h>

#include "ggsw_ct_params.h"

/** @brief opaque type that are included in */
typedef struct module_info_t MODULE;

/** @brief opaque type that represents a vector of znx in DFT space */
typedef struct vec_znx_dft_t VEC_ZNX_DFT;

/** @brief opaque type that represents a prepared matrix */
typedef struct vmp_pmat_t VMP_PMAT;


/** @brief opaque type that represents a bivariate polynomial in DFT space 
 *  @note `Length = l*N`
*/
typedef double PolyBivDFT;

/** @brief opaque type that represents a bivariate polynomial vector in DFT space 
 *  @note `Length = (k+1) * l * N` for a bivGLWE
*/
typedef double VecBivDFT;

/** @brief opaque type that represents a bivariate polynomial matrix in DFT space
 *  @note `Length = n_limbs_tilde * (k+1) * l * N` for a bivGGSW
*/
typedef double MatBivDFT;

/** @brief opaque type that represents a univariate polynomial in DFT space 
 *  @note `Length = N for a polynomial in ZnX`
*/
typedef double PolyUnivDFT;

/** @brief opaque type that represents a univariate polynomial in DFT space 
 *  @note `Length = k * N for a GLWE`
*/
typedef double VecUnivDFT;


/** @brief opaque type that represents a bivariate polynomial 
 *  @note `Length = l * N for bivariate polynomial in ZnXY`
*/
typedef double PolyBiv;

/** @brief opaque type that represents a bivariate polynomial vector 
 *  @note `Length = l * k * N for a bivGLWE`
*/
typedef double VecBiv;

/** @brief opaque type that represents a bivariate polynomial matrix 
 *  @note `Length = n_limbs_tilde * l * k * N
*/
typedef double MatBiv;

/** @brief opaque type that represents a univariate polynomial 
 *  @note `Length = N for a polynomial in ZnX`
*/
typedef double PolyUniv;


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
  MatBiv* ct; //represent a matrix of size n_limbs x n_limbs_tilde with coefficients that are in ZnX
} GGSWCiphertext;

typedef struct ggsw_ciphertext_dft {
  GGSWCtParams* params;
  MatBivDFT* ct;
} GGSWPreparedCt;

#endif  // GGSW_CIPHERTEXT_H
