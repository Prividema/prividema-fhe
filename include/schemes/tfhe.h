#ifndef TFHE_H
#define TFHE_H

#include "ggsw_ciphertext.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"

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
 *                      Othwerise it is assumed that c0 and c1 are close enough to normal form that the external
 *                      product will behave correctly.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int tfhe_cmux_unprepared(const PvdaBackend* module, GLWECiphertext* res, const GLWECiphertext* c0,
                         const GLWECiphertext* c1, const GGSWCiphertext* sel, int normalize_sub);

/**
 * @brief Performs a CMux (secret inputs and selector multiplex) operation
 *
 * @param module The backend module
 * @param res The output GLWE, which is equivalent to c0 if sel encrypts 0, to c1 if sel encrypts 1,
 *            and undefined otherwise
 * @param c0 Input for the multiplexer selected with signal 0
 * @param c1 Input for the multiplexer selected with signal 1
 * @param sel The selector prepared GGSW ciphertext
 * @param normalize_sub Boolean, if true we normalize after subtracting c0 and c1.
 *                      Othwerise it is assumed that c0 and c1 are close enough to normal form that the external
 *                      product will behave correctly.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int tfhe_cmux(const PvdaBackend* module, GLWECiphertext* res, const GLWECiphertext* c0, const GLWECiphertext* c1,
              const GGSWCiphertextPrep* sel, int normalize_sub);

/**
 * @brief Select one GLWE ciphertext out of many using an encrypted selection signal, using a Mux tree
 *
 * Implements a CMux tree to select between inp_cols GLWEs provided in src.
 * It takes as selection signals GGSWCiphertextPreps that encode either 0 or 1.
 * selectors[0] encodes the LSB of the column number that is to be selected,
 * selectros[1] the second least significant bit of the column number, etc.
 * Column numbers start from 0.
 * sel_size must be at least ceil(log2(inp_cols))
 *
 * It can optionally deallocate the input ciphertexts for reduced memory consumption
 *
 * @param module The backend module
 * @param res The GLWE result
 * @param src A vector of GLWECiphertexts to be selected from
 * @param inp_cols The size of the vector, what is, how many GLWEs to select from
 * @param selectors The selection signals
 * @param sel_size The number of provided selection signals
 * @param delete_src Whether to delete the input GLWEs, 1 to enable it and 0 to preserve the input
 *
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int tfhe_cmux_tree(const PvdaBackend* module, GLWECiphertext* res, const GLWECiphertext** src, int inp_cols,
                   const GGSWCiphertextPrep** selectors, int sel_size, int delete_src);

void tfhe_blindrotate(PvdaBackend* module, GLWECiphertext* res, const GLWECiphertext*);

#endif
