#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

#include <stdint.h>

#include "ggsw_encrypt_params.h"
#include "halfggsw_ciphertext.h"

typedef struct ggsw_ciphertext {
  GGSWEncryptParams* params;
  int64_t* ct; //represent a matrix of size n_limbs x n_limbs_tilde with coefficients that are in ZnX
    //  values;  // BASE2K_INT_POL = ZNX ? Vector of RLWE ? ? Multiple halfggsw
  // values[x][y] takes the x-th line, y-th column of the ciphertext
  // typically, there are (k+1).l lines and k columns
} GGSWCiphertext;


#endif  // GGSW_CIPHERTEXT_H
