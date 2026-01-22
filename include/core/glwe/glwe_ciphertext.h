#ifndef GLWE_CIPHERTEXT_H
#define GLWE_CIPHERTEXT_H

#include <stdint.h>
#include "glwe_ct_params.h"
#include "bivariate_polynomial.h"
#include "spqlios_alias.h"

//! GLWE PART (begin)

typedef struct glwe_ciphertext {
  GLWECtParams* params;   // GLWE parameters
  VecBiv* vec;            // Represents a vector of size (k + 1) * l with coefficients that are in Zn[X]
} GLWECiphertext;

/**
 * @brief Return the number of coefficient in a bivariate GLWE ciphertext.
 * 
 * @param params The GLWE parameters.
 * @return int64_t 
 */
int64_t glwe_coef_number(GLWECtParams* params);

/**
 * @brief Creates a bivGLWE, filled with 0.
 * 
 * @param params The GLWE parameters.
 * @return GLWECiphertext*
 */
GLWECiphertext* new_glwe(GLWECtParams* params);
void delete_glwe(GLWECiphertext* ct);
void add_glwe(GLWECiphertext* res, GLWECiphertext* ct1, GLWECiphertext* ct2);
void const_mult_glwe(GLWECiphertext* res, PolyUniv* u, GLWECiphertext* ct);


//! GLWE IN DFT PART (begin)

typedef struct glwe_ciphertext_dft {
  GLWECtParams* params;     // GLWE parameters
  VecBivDFT* pvec;          // Prepared vector
} GLWEPreparedCt;

/**
 * @brief The number of coefficient in a bivariate GLWE ciphertext in DFT space.
 * 
 * @param params The GLWE parameters.
 * @return int64_t 
 * 
 * @note The number of independent coefficients of a polynomial in DFT space is half the number of coefficients in Zn[X], 
 * due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
int64_t glwe_coef_number_dft(GLWECtParams* params);
GLWEPreparedCt* new_glwe_dft(GLWECtParams* params);
void delete_glwe_dft(GLWEPreparedCt* ct);
void add_glwe_dft(GLWEPreparedCt* res, GLWEPreparedCt* ct1, GLWEPreparedCt* ct2);
void const_mult_glwe_dft(GLWEPreparedCt* res, PolyUniv* u, GLWEPreparedCt* ct);

//! COMMON PART (begin)

/**
 * @brief Return the size of a bivGLWE ciphertext, in DFT space & out of DFT space.
 * 
 * @param params The GLWE parameters.
 * @return int64_t 
 * 
 * @note The size of a bivGLWE ciphertext is the same in and out of DFT space.
 */
int64_t glwe_size(GLWECtParams* params);

/**
 * @brief The number of bytes needed to store a bivGLWE ciphertext.
 * 
 * @param params The GLWE parameters.
 * @return int64_t 
 * 
 * @note The number of bytes needed to store a bivGLWE ciphertext, is the same in and out of DFT space. 
 */
int64_t glwe_bytes(GLWECtParams* params);

/**
 * @brief Compute the polynomial product of c and d, component-wise in DFT space.
 * 
 * @param module The module stocking the degree N.
 * @param res_dft The result in DFT space.
 * @param res_size The result's size.
 * @param c_dft The left-hand side polynomial in DFT space .
 * @param c_size The left-hand size of c_dft.
 * @param d_dft The right-hand side polynomial in DFT space.
 * @param d_size The right-hand size of c_dft.
 * 
 * @note `res_dft = ( DFT(c_0) * DFT(d_0) , ... , DFT(c_smin) * DFT(d_smin) , 0's)`. There are enough 0's to match the size of res_dft.
 */
void mult_vec_znx_dft(const MODULE* module, 
              double* res_dft, int64_t res_size,
              double* c_dft, int64_t c_size,  
              double* d_dft, int64_t d_size
);


#endif  // GLWE_CIPHERTEXT_H