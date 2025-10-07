#ifndef HALFGGSW_CIPHERTEXT_H
#define HALFGGSW_CIPHERTEXT_H

#include "ggsw_ciphertext.h"
#include "glwe_ciphertext.h"

typedef struct halfggsw_ciphertext{
    GGSW_ENCRYPT_PARAMS* params;
    GLWE_CIPHERTEXT* values; // vector of RLWE
} HALFGGSW_CIPHERTEXT;

#endif // HALFGGSW_CIPHERTEXT_H