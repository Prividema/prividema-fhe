#ifndef TFHE_H
#define TFHE_H

#include "ggsw_arithmetic.h"
#include "glwe_arithmetic.h"

/**
 * @file tfhe.h
 * @brief TFHE-related operations
 *
 */

/**
 * @brief Performs a CMux (secret inputs and selector multiplex) operation
 *
 * @param module The backend module
 * @param res The output GLWE, which is equivalent to c0 if sel encrypts 0, to c1 if sel encrypts 1,
 *            and undefined otherwise
 * @param c0 Input for the multiplexer selected with signal 0
 * @param c1 Input for the multiplexer selected with signal 1
 * @param sel The selector GGSW ciphertext
 * @param normalize_sub Boolean, if true we normalize after subtracting c0 and c1.
 *                      Othwerise we hope that c0 and c1 are close enough to normal form that the external
 *                      product will behave correctly.
 *
 */
int tfhe_cmux_unprepared(MODULE* module, GLWECiphertext* res, const GLWECiphertext* c0, const GLWECiphertext* c1,
                         const GGSWCiphertext* sel, int normalize_sub);

int tfhe_cmux(MODULE* module, GLWECiphertext* res, const GLWECiphertext* c0, const GLWECiphertext* c1,
              const GGSWCiphertextPrep* sel, int normalize_sub);

void tfhe_blindrotate(MODULE* module, GLWECiphertext* res, const GLWECiphertext*);

#endif
