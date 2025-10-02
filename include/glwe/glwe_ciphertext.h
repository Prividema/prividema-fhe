#ifndef GLWE_CIPHERTEXT_H
#define GLWE_CIPHERTEXT_H

#include <stdint.h>

typedef struct glwe_ciphertext {
    uint64_t N;
    uint64_t k;
    uint64_t kappa; //base_2k
    uint64_t l_a; //size of the decomposition of a
    uint64_t l_b; //size of the decomposition of b
    int64_t* b; 
    int64_t* a; // vector of k elements
    void* data; // memory of the struct
} GLWE_CIPHERTEXT;

#endif // GLWE_CIPHERTEXT_H