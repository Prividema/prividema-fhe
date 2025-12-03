#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

#include <stdint.h>

#include "ggsw_ct_params.h"

typedef struct ggsw_ciphertext {
  GGSWCtParams* params;

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

  int64_t* ct; //represent a matrix of size n_limbs x n_limbs_tilde with coefficients that are in ZnX
} GGSWCiphertext;


#endif  // GGSW_CIPHERTEXT_H
