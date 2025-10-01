#ifndef HALFGGSW_CIPHERTEXT_H
#define HALFGGSW_CIPHERTEXT_H

#include <stdint.h>
#include "ggsw_ciphertext.h"

typedef struct halfgsw_ciphertext{
    GSW_ENCRYPT_PARAMS* params;
    BASE2K_INT_POL** values; // BASE2K_INT_POL = ZNX ? Maybe vector of RLWE ?
    // values[x][y] takes the x-th line, y-th column of the ciphertext
    // typically, there are l lines and k columns
} HALFGSW_CIPHERTEXT;

#endif