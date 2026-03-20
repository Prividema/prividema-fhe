#ifndef TFHE_H
#define TFHE_H

#include "ggsw.h"
#include "glwe.h"

/*
    Implements a homomorphic selection (Mux) operation,
    where GGSW ciphertext sel is used to choose from c0 or c1

    TODO: finish
*/
void tfhe_cmux(MODULE* module, GLWECiphertext* res, const GLWECiphertext* c0, const GLWECiphertext* c1,
               const GGSWCiphertext* sel);

void tfhe_blindrotate(MODULE* module, GLWECiphertext* res, const GLWECiphertext*);

#endif
