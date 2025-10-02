#ifndef GLWE_KEY_H
#define GLWE_KEY_H

#include "glwe_ciphertext.h"

// TODO : To be Defined (spqlios).
typedef struct core {} CORE;
typedef struct svp_ppol {} SVP_PPOL;

typedef struct glwe_secret_key {
    uint64_t n; // n= N*k (n is constant: N and k can change when we change the module N) 
    int64_t* s; 
    void* data; 
} GLWE_SECRET_KEY;

typedef struct glwe_prep_secret_key {
    uint64_t N;
    uint64_t k;
    SVP_PPOL** s; // vec of size k, each element is prepared vec
    void* data;
} GLWE_PREP_SECRET_KEY;

typedef struct glwe_public_key {
    uint64_t N;
    uint64_t k;
    uint64_t Y;
    GLWE_CIPHERTEXT** pk; // vector of Y element (A1,...,Ak, B)
    void* data; 
} GLWE_PUBLIC_KEY;

//generation secret key
void glwe_secret_key_gen (
                        const CORE* core,
                        uint64_t lambda, 
                        GLWE_SECRET_KEY* s
);

//generation public key
void glwe_public_key_gen (
                        const CORE* core,
                        uint64_t lambda, 
                        GLWE_SECRET_KEY* s,
                        GLWE_PUBLIC_KEY* pk
);


#endif // GLWE_KEY_H