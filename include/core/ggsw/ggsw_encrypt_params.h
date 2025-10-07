#ifndef GGSW_ENCRYPT_PARAMS_H
#define GGSW_ENCRYPT_PARAMS_H

#include <stdint.h>
#include "glwe_encrypt_params.h"

typedef struct ggsw_encrypt_params{
    GLWE_ENCRYPT_PARAMS* params;
    int64_t l; // Decomposition size for the gadget
    int64_t B; // Decomposition basis for the gadget
} GGSW_ENCRYPT_PARAMS;

#endif // GGSW_ENCRYPT_PARAMS_H