#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

#include <stdint.h>

#include "ggsw_encrypt_params.h"
#include "glwe.h"
#include "vec_znx_arithmetic_private.h"
#include "halfggsw_encrypt_params.h"

/** @brief Defintion of HalfBivGGSW ciphertexts */

uint64_t index_halfggsw_vmp(uint64_t k, uint64_t DegreeInY, uint64_t N, 
  uint64_t i, uint64_t j, uint64_t jj, uint64_t jjj
){
  return (i*(k+1)*DegreeInY + j*DegreeInY + jj)*N + jjj;
}
typedef typeof(index_halfggsw_vmp) INDEX_HALFGGSW;

typedef VMP_PMAT VecBivRLWE;

typedef struct halfggsw_ciphertext {

  INDEX_HALFGGSW* index;
  HalfGGSWEncryptParams* params; // Parameters for the encryption and external/internal product
  VecBivRLWE* values;  // Prepared vector of poly representing a vector of bivRLWE
  
} HalfGGSWCiphertext;

typedef struct ggsw_ciphertext {

  GGSWEncryptParams* params; // Parameters for the encryption and external/internal product
  HalfGGSWCiphertext** halfggswVec; // vector of k HalfGGSW Ciphertexts

} GGSWCiphertext;

HalfGGSWCiphertext* new_halfggsw(Core* core, uint64_t nrows, uint64_t DegreeInY); 

GGSWCiphertext* new_ggsw(Core* core, uint64_t K, uint64_t K_tilde, uint64_t L, uint64_t l_tilde);
    
#endif  // GGSW_CIPHERTEXT_H
