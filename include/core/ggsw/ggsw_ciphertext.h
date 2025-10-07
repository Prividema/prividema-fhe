#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

#include <stdint.h>
#include "ggsw_encrypt_params.h"
#include "halfggsw_ciphertext.h"

typedef struct ggsw_ciphertext{
    GGSW_ENCRYPT_PARAMS* params;
    HALFGGSW_CIPHERTEXT* values; // BASE2K_INT_POL = ZNX ? Vector of RLWE ? ? Multiple halfggsw
    // values[x][y] takes the x-th line, y-th column of the ciphertext
    // typically, there are (k+1).l lines and k columns
} GGSW_CIPHERTEXT;

#endif // GGSW_CIPHERTEXT_H