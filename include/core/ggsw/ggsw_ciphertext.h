#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

/**
 * @file ggsw_ciphertext.h
 * @brief GGSW ciphertext
 * 
 * In this header file, we define the structure representing bivariate GGSW ciphertext 
 * and bivariate GGSW in DFT space. Moreover, we define alias for double and int64_t. 
 * Thus, we properly manipulate mathematical structure (polynomial, vector or matrix).
 * The types below are provided solely for clarity, nothing is changed in spqlios.
 */

#include <stdint.h>
#include <stdlib.h>
#include "ggsw_ct_params.h"
#include "structure_alias.h"
#include "partialggsw_ciphertext.h"


//! GGSW PART (begin)

typedef struct ggsw_ciphertext {
  GGSWCtParams* params;  // GGSW parameters
  MatBiv* mat;           // Represent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in ZnX
} GGSWCiphertext;

int new_ggsw(GGSWCiphertext* res, GGSWCtParams* params, MatBiv* ct);
void delete_ggsw(GGSWCiphertext* ct);
VecBiv* ggsw_Sj_Yi(GGSWCiphertext* ct, int64_t i, int64_t j);
void add_ggsw(GGSWCiphertext* res, GGSWCiphertext* ct1, GGSWCiphertext* ct2);
void const_mult_ggsw(GGSWCiphertext* res, GGSWCiphertext* ct, PolyUniv* u);


//! GGSW IN DFT PART
typedef struct ggsw_ciphertext_dft {
  GGSWCtParams* params;  // GGSW parameters
  MatBivDFT* pmat;       // Represent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in ZnX
} GGSWCiphertextDFT;

int new_ggsw_dft(GGSWCiphertextDFT* res, GGSWCtParams* params, MatBivDFT* ct);
void delete_ggsw_dft(GGSWCiphertextDFT* ct);
VecBivDFT* ggsw_Sj_Yi_dft(GGSWCiphertextDFT* ct, int64_t i, int64_t j);
void add_ggsw_dft(GGSWCiphertext* res, GGSWCiphertext* ct1, GGSWCiphertext* ct2);
void const_mult_ggsw_dft(GGSWCiphertextDFT* res_dft, GGSWCiphertextDFT* ct_dft, PolyUniv* u);

//! COMMON PART (begin)

int64_t ggsw_size(GGSWCtParams* params);
int64_t ggsw_bytes(GGSWCtParams* params);


#endif  // GGSW_CIPHERTEXT_H
