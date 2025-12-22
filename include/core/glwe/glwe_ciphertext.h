#ifndef GLWE_CIPHERTEXT_H
#define GLWE_CIPHERTEXT_H

#include <stdint.h>
#include "glwe_ct_params.h"
#include "structure_alias.h"


typedef struct glwe_ciphertext {
  GLWECtParams* params;
} GLWECiphertext;

//! Basics functions (begin)

int new_glwe(GLWECiphertext* res, GLWECtParams* params);
void delete_glwe(GLWECiphertext* ct);
void add_glwe(GLWECiphertext* res, GLWECiphertext* ct1, GLWECiphertext* ct2);
void const_mult_glwe(GLWECiphertext* res, PolyUniv* u, GLWECiphertext* ct);

typedef struct glwe_ciphertext_dft {
  GLWECtParams* params;
} GLWEPreparedCt;

//! Basics functions (begin)

int new_glwe_dft(GLWEPreparedCt* res, GLWECtParams* params);
void delete_glwe_dft(GLWEPreparedCt* ct);
void add_glwe_dft(GLWEPreparedCt* res, GLWEPreparedCt* ct1, GLWEPreparedCt* ct2);
void const_mult_glwe_dft(GLWEPreparedCt* res, PolyUniv* u, GLWEPreparedCt* ct);

#endif  // GLWE_CIPHERTEXT_H
