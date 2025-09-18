#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

#include <cstdint>
#include "ggsw_encrypt_params.h"
//#include something for BASE2K_INT_POL

/* /!\ Should it be represented as two halfggsw ? */
typedef struct gsw_ciphertext{
    GSW_ENCRYPT_PARAMS* params;
    BASE2K_INT_POL** values; // BASE2K_INT_POL = ZNX ?
    // values[x][y] takes the x-th line, y-th column of the ciphertext
    // typically, there are (k+1).l lines and k columns
} GSW_CIPHERTEXT;

#endif