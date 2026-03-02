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

typedef struct ggsw_ciphertext
{
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
 * @brief Creates a bivGGSW, if ct = NULL it is filled with 0, otherwise the bivGGSW ciphertext is a copy a ct.
 *
 * @param params The GGSW parameters.
 * @param ct The GGSW matrix.
 *
 * @retval - `NULL` if malloc failed inside the function.
 * @retval - Allocated GGSW ciphertext, filled with 0 or ct, if ct is non-null. 
 * 
 */
GGSWCiphertext* new_ggsw(GGSWCtParams* params, MatBiv* ct);

/**
 * @brief Delete a bivGGSW ciphertext.
 *
 * @param ct The bivGGSW ciphertext.
 */
void delete_ggsw(GGSWCiphertext* ct);

/**
 * @brief Return the pointer to biGLWE(-m * sk_j / (2^kappa_tilde)^(i+1)).
 *
 * @param params_ggsw The GGSW parameters.
 * @param ct_mat The GGSW ciphertext's matrix.
 * @param j The j-th component of Sk.
 * @param i The i in the phase = -m * sk_j / (2^kappa_tilde)^(i+1).
 *
 * @return VecBiv*
 */
VecBiv* ggsw_Sj_Yti(GGSWCtParams* params_ggsw, MatBiv* ct_mat, int64_t j, int64_t i);

/**
 * @brief Normalize a GGSW ciphertext.
 *
 * @param res The result normalized GGSW ciphertext.
 * @param ct The input GGSW ciphertext.
 */
void normalize_ggsw(MODULE* module, GGSWCiphertext* res, GGSWCiphertext* ct);

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
 * @param u_dft The Zn[X] polynomial in DFT space.
 */
int const_mult_ggsw(MODULE* module, GGSWCiphertext* res, GGSWCiphertext* ct, PolyUnivDFT* u, int do_normalization);

//! GGSW IN DFT PART
// The same functions but the operation are done in DFT space. The parameters are the same but the ciphertext are in DFT
// space.

typedef struct ggsw_ciphertext_dft
{
	GGSWCtParams* params;  // GGSW parameters
	MatBivDFT* mat;        // Represent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in Zn[X]
} GGSWCiphertextDFT;

/**
 * @brief The number of coefficient in a bivariate GGSW ciphertext in DFT space.
 *
 * @param params The GGSW parameters.
 * @return int64_t
 *
 * @note The number of independent coefficients of a polynomial in DFT space is half the number of coefficients in
 * Zn[X], due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t ggsw_coef_number_dft(GGSWCtParams* params);

/**
 * @brief Creates a bivGGSW in DFT space, filled with 0.
 *
 * @param params The GGSW parameters.
 * @param mat_dft The GGSW matrix in DFT space.
 */
GGSWCiphertextDFT* new_ggsw_dft(GGSWCtParams* params, MatBivDFT* ct);

/**
 * @brief Delete a bivGGSW in DFT space.
 *
 * @param ct The bivGGSW in DFT space.
 */
void delete_ggsw_dft(GGSWCiphertextDFT* ct);

/**
 * @brief Return the pointer to biGLWE(DFT(-m * sk_j / 2^kappa*(i+1))) in DFT space.
 *
 * @param params_ggsw The GGSW parameters.
 * @param ct_mat The GGSW ciphertext's matrix.
 * @param j The j-th component of Sk.
 * @param i The i in the phase = -m * sk_j / (2^kappa_tilde)^(i+1).
 *
 * @return VecBivDFT*
 */
VecBivDFT* ggsw_Sj_Yti_dft(GGSWCtParams* params_ggsw, MatBivDFT* ct_mat, int64_t j, int64_t i);

/**
 * @brief Normalize a GGSW ciphertext in DFT space.
 *
 * @param res_dft The result normalized GGSW ciphertext in DFT space.
 * @param ct_dft The input GGSW ciphertext in DFT space.
 */
int normalize_ggsw_dft(MODULE* module, GGSWCiphertextDFT* res_dft, GGSWCiphertextDFT* ct_dft);

/**
 * @brief Adds two GGSW ciphertexts in DFT space with same params and put result in res.
 *
 * @param res The result GGSW ciphertext in DFT space.
 * @param ct1 The left-hand side GGSW ciphertext in DFT space.
 * @param ct2 The right-hand side GGSW ciphertext in DFT space.
 */
void add_ggsw_dft(GGSWCiphertextDFT* res_dft, GGSWCiphertextDFT* ct1_dft, GGSWCiphertextDFT* ct2_dft);

/**
 * @brief  Multiply a GGSW ciphertext by a constant in Zn[X]
 *
 * @param res_dft The result GGSW ciphertext.
 * @param ct_dft The GGSW ciphertext.
 * @param u_dft The Zn[X] polynomial in DFT space, with coefficient in [-2^(kappa-1), 2^(kappa-1)]
 */
int const_mult_ggsw_dft(MODULE* module, GGSWCiphertextDFT* res_dft, GGSWCiphertextDFT* ct_dft, PolyUnivDFT* u_dft,
                        int do_normalization);

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