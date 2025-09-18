#ifndef GGSW_PUBLIC_KEY_H
#define GGSW_PUBLIC_KEY_H

#include <cstdint>
#include "ggsw_ciphertext.h"

typedef struct gsw_public_key{
    GSW_CIPHERTEXT* pk; // Public key is multiple encryptions of 0
} GSW_PUBLIC_KEY;

#endif