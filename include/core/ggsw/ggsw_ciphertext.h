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
#include "maths_structures.h"
#include "partialggsw_ciphertext.h"


//! GGSW PART (begin)

typedef struct ggsw_ciphertext {
    GGSWCtParams* params;  // GGSW parameters
    MatBiv* mat;           // Represent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in Zn[X]
} GGSWCiphertext;


/**
 * @brief The number of coefficients in a bivariate GGSW ciphertext.
 * 
 * @param params The GGSW parameters.
 * @return int64_t 
 */
uint64_t ggsw_coef_number(GGSWCtParams* params);

/**
 * @brief Creates a bivGGSW, filled with 0.
 * 
 * @param params The GGSW parameters.
 * @param mat The GGSW matrix.
 * 
 * @return GGSWCiphertext*
 */
GGSWCiphertext* new_ggsw(GGSWCtParams* params, MatBiv* ct);

/**
 * @brief Delete a bivGGSW ciphertext.
 * 
 * @param ct The bivGGSW ciphertext.
 */
void delete_ggsw(GGSWCiphertext* ct);

/**
 * @brief Return the pointer to biGLWE(-m * sk_j / (2^kappa_tilde)^i).
 * 
 * @param ct A GGSW ciphertext.
 * @param i The degree in Y of the phase = -m * sk_j / (2^kappa_tilde)^i.
 * @param j The j-th component of Sk.
 * 
 * @return VecBiv*
 */
VecBiv* ggsw_Sj_Yti(GGSWCiphertext* ct, int64_t i, int64_t j);

/**
 * @brief Normalize a GGSW ciphertext.
 * 
 * @param res The result normalized GGSW ciphertext.
 * @param ct The input GGSW ciphertext.
 */
void normalize_ggsw(GGSWCiphertext* res, GGSWCiphertext* ct);

/**
 * @brief Adds two GGSW ciphertexts with same params and put result in res
 * 
 * @param res The result GGSW ciphertext. 
 * @param ct1 The left-hand side GGSW ciphertext.
 * @param ct2 The right-hand side GGSW ciphertext.
 */
void add_ggsw(GGSWCiphertext* res, GGSWCiphertext* ct1, GGSWCiphertext* ct2);

/**
 * @brief  Multiply a GGSW ciphertext by a constant in Zn[X]
 * 
 * @param res The result GGSW ciphertext.
 * @param ct The GGSW ciphertext.
 * @param u The polynomial in Zn[X], with coefficient in [-2^(kappa-1), 2^(kappa-1)]
 */
void const_mult_ggsw(GGSWCiphertext* res, GGSWCiphertext* ct, PolyUniv* u, int do_normalization);


//! GGSW IN DFT PART
// The same functions but the operation are done in DFT space. The parameters are the same but the ciphertext are in DFT space.

typedef struct ggsw_ciphertext_dft {
    GGSWCtParams* params;  // GGSW parameters
    MatBivDFT* pmat;       // Represent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in Zn[X]
} GGSWCiphertextDFT;

/**
 * @brief The number of coefficient in a bivariate GGSW ciphertext in DFT space.
 * 
 * @param params The GGSW parameters.
 * @return int64_t 
 * 
 * @note The number of independent coefficients of a polynomial in DFT space is half the number of coefficients in Zn[X], 
 * due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t ggsw_coef_number_dft(GGSWCtParams* params);

/**
 * @brief Creates a bivGGSW in DFT space, filled with 0.
 * 
 * @param params The GGSW parameters.
 * @param mat_dft The GGSW matrix in DFT space.
 */
GGSWCiphertextDFT* new_ggsw_dft( GGSWCtParams* params, MatBivDFT* ct);

/**
 * @brief Delete a bivGGSW in DFT space.
 * 
 * @param ct The bivGGSW in DFT space.
 */
void delete_ggsw_dft(GGSWCiphertextDFT* ct);

/**
 * @brief Return the pointer to biGLWE(DFT(-m * sk_j) * Y^i) in DFT space.
 * 
 * @param ct_dft A GGSW ciphertext in DFT space.
 * @param i The degree in Y of the phase = -m * sk_j * Y^i.
 * @param j The j-th component of Sk.
 * 
 * @return VecBivDFT*
 */
VecBivDFT* ggsw_Sj_Yti_dft(GGSWCiphertextDFT* ct, int64_t i, int64_t j);

/**
 * TODO : Implement it
 */
void add_ggsw_dft(GGSWCiphertext* res, GGSWCiphertext* ct1, GGSWCiphertext* ct2);

/**
 * @brief  Multiply a GGSW ciphertext by a constant in Zn[X]
 * 
 * @param res_dft The result GGSW ciphertext.
 * @param ct_dft The GGSW ciphertext.
 * @param u The polynomial in Zn[X], with coefficient in [-2^(kappa-1), 2^(kappa-1)]
 */
void const_mult_ggsw_dft(GGSWCiphertextDFT* res_dft, GGSWCiphertextDFT* ct_dft, PolyUniv* u);

//! COMMON PART (begin)

/**
 * @brief Return the size of a bivGGSW ciphertext, in DFT space & out of DFT space.
 * 
 * @param params The GGSW parameters.
 * @return int64_t 
 * 
 * @note The size of a bivGGSW ciphertext is the same in and out of DFT space.
 */
uint64_t ggsw_size(GGSWCtParams* params);

/**
 * @brief The number of bytes needed to store a bivGGSW ciphertext.
 * 
 * @param params The GGSW parameters.
 * @return int64_t 
 * 
 * @note The number of bytes needed to store a bivGGSW ciphertext, is the same in and out of DFT space. 
 */
uint64_t ggsw_bytes(GGSWCtParams* params);


#endif  // GGSW_CIPHERTEXT_H