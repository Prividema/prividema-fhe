#ifndef GGSW_PUBLIC_KEY_H
#define GGSW_PUBLIC_KEY_H

#include <stdint.h>
#include "ggsw_ciphertext.h"

typedef struct gsw_public_key{
    GSW_CIPHERTEXT* pk; // Public key is multiple encryptions of 0
    // Another thing to count size of key ?
} GSW_PUBLIC_KEY;

#endif // GGSW_PUBLIC_KEY_H