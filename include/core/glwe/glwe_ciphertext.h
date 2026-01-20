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


int64_t glwe_coef_number(GLWECtParams* params);
GLWECiphertext* new_glwe(GLWECtParams* params);
void delete_glwe(GLWECiphertext* ct);
void add_glwe(GLWECiphertext* res, GLWECiphertext* ct1, GLWECiphertext* ct2);
void const_mult_glwe(GLWECiphertext* res, PolyUniv* u, GLWECiphertext* ct);


//! GLWE IN DFT PART (begin)

typedef struct glwe_ciphertext_dft {
  GLWECtParams* params;     // GLWE parameters
  VecBivDFT* pvec;          // Prepared vector
} GLWEPreparedCt;


int64_t glwe_coef_number_dft(GLWECtParams* params);
int new_glwe_dft(GLWEPreparedCt* res, GLWECtParams* params);
void delete_glwe_dft(GLWEPreparedCt* ct);
void add_glwe_dft(GLWEPreparedCt* res, GLWEPreparedCt* ct1, GLWEPreparedCt* ct2);
void const_mult_glwe_dft(GLWEPreparedCt* res, PolyUniv* u, GLWEPreparedCt* ct);

//! COMMON PART (begin)

int64_t glwe_size(GLWECtParams* params);
int64_t glwe_bytes(GLWECtParams* params);

void vec_znx_dft_mult(const MODULE* module, 
              double* res_dft, int64_t res_size,
              double* c_dft, int64_t c_size,  
              double* d_dft, int64_t d_size
);


#endif  // GLWE_CIPHERTEXT_H