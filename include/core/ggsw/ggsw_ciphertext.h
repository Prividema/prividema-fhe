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


typedef struct ggsw_ciphertext {
  GGSWCtParams* params;
  MatBiv* ct; //represent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in ZnX
} GGSWCiphertext;

int new_ggsw(GGSWCiphertext* res, GGSWCtParams* params, MatBiv* ct);
void delete_ggsw(GGSWCiphertext* ct);
void add_ggsw(GGSWCiphertext* res, GGSWCiphertext* ct1, GGSWCiphertext* ct2);
void const_mult_ggsw();


typedef struct ggsw_ciphertext_dft {
  GGSWCtParams* params;
  MatBivDFT* ct;
} GGSWPreparedCt;

int new_ggsw_prepared(GGSWPreparedCt* res, GGSWCtParams* params, MatBivDFT* ct);
void delete_ggsw_prepared(GGSWPreparedCt* ct);
void add_ggsw_prepared(GGSWCiphertext* res, GGSWCiphertext* ct1, GGSWCiphertext* ct2);

#endif  // GGSW_CIPHERTEXT_H
