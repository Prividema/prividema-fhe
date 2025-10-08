#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

#include <stdint.h>

#include "ggsw_encrypt_params.h"
#include "halfggsw_ciphertext.h"

typedef struct ggsw_ciphertext {
  GGSWEncryptParams* params;
  HalfGGSWCiphertext*
      values;  // BASE2K_INT_POL = ZNX ? Vector of RLWE ? ? Multiple halfggsw
  // values[x][y] takes the x-th line, y-th column of the ciphertext
  // typically, there are (k+1).l lines and k columns
} GGSWCiphertext;

#endif  // GGSW_CIPHERTEXT_H
