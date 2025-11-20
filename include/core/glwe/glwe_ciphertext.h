#ifndef GLWE_CIPHERTEXT_H
#define GLWE_CIPHERTEXT_H

#include <stdint.h>
#include "glwe_ct_params.h"

// TODO : To be Defined (spqlios).
typedef struct b {
} Base2kIntegerPoly;

typedef struct glwe_ciphertext {
  GLWECtParams* params;
  int64_t* ct;  //  corresponds to (a_0, ..., a_{k-1}, b=a_k)
  
  // l_a = floor((n_limbs+1)/(k+1))
  // l_b= n_limbs- k*l_a
  // L (logalpha) = l_b*kappa //precision
  // ct[p]= a_i[j], i= p%(k+1) in [0, l_a[ and j= floor(p/k+1) in [O,k[
  // a0(Y) = a0[O] + a0[1]Y + .... + a0[la-1]Y^{la-1}
  //
  // ak(Y) = ak[0]+ ak[1]Y + ... + ak[lb-1]Y^{lb-1}
} GLWECiphertext;

#endif  // GLWE_CIPHERTEXT_H
