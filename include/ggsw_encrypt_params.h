#ifndef GGSW_ENCRYPT_PARAMS_H
#define GGSW_ENCRYPT_PARAMS_H

#include <cstdint>

typedef struct gsw_encrypt_params{
    int64_t N; // Polynomial degree
    int64_t k; // Number of a terms
    int64_t l; // Decomposition size for the gadget
    int64_t B; // Decomposition basis for the gadget
    int64_t base; // Base used for the base-2^k representation
    int64_t prec; // Decomposition size for base-2^k representation
} GSW_ENCRYPT_PARAMS;

#endif