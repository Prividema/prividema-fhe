#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

/**
 * @file ggsw_ciphertext.h
 * @brief GGSW ciphertext
 *
 * In this header file, we define the structure represultenting bivariate GGSW ciphertext
 * and bivariate GGSW in the DFT domain. Moreover, we define alias for double and int64_t.
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
	const GGSWCtParams* params;  // GGSW parameters
	MatBiv* mat;  // Represultent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in Zn[X]
} GGSWCiphertext;

/**
 * @brief The number of coefficients in a bivariate GGSW ciphertext.
 *
 * @param params_ggsw The GGSW parameters.
 * @return int64_t
 */
uint64_t ggsw_coef_number(const GGSWCtParams* params_ggsw);

/**
 * @brief Creates a bivGGSW, it is filled with 0.
 *
 * @param params_ggsw The GGSW parameters.
 *
 * @retval - `NULL` if malloc failed inside the function.
 * @retval - Allocated GGSW ciphertext, filled with 0.
 *
 */
GGSWCiphertext* new_ggsw(const GGSWCtParams* params_ggsw);

/**
 * @brief Delete a bivGGSW ciphertext.
 *
 * @param ggsw The bivGGSW ciphertext.
 */
void delete_ggsw(GGSWCiphertext* ggsw);

/**
 * @brief Return the pointer to biGLWE(-m * sk_j / (2^kappa_tilde)^(i+1)).
 *
 * @param params_ggsw The GGSW parameters.
 * @param ggsw_mat The GGSW ciphertext's matrix.
 * @param j The j-th component of Sk.
 * @param i The i in the phase = -m * sk_j / (2^kappa_tilde)^(i+1).
 *
 * @return VecBiv*
 */
VecBiv* ggsw_Sj_Yti(const GGSWCtParams* params_ggsw, MatBiv* ggsw_mat, int64_t j, int64_t i);

/**
 * @brief Normalize a GGSW ciphertext.
 *
 * @param module Additionnal information for backend.
 * @param result The result normalized GGSW ciphertext.
 * @param ggsw The input GGSW ciphertext.
 */
int normalize_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw);

/**
 * @brief Adds two GGSW ciphertexts with same params and put result in result
 *
 * @param result The result GGSW ciphertext.
 * @param ggsw_lhs The left-hand side GGSW ciphertext.
 * @param ggsw_rhs The right-hand side GGSW ciphertext.
 */
void add_ggsw(GGSWCiphertext* result, const GGSWCiphertext* ggsw_lhs, const GGSWCiphertext* ggsw_rhs);

/**
 * @brief Multiply a GGSW ciphertext by a constant in Zn[X]
 *
 * @param module Additionnal information for backend.
 * @param result The result GGSW ciphertext.
 * @param ggsw The GGSW ciphertext.
 * @param u_dft The Zn[X] polynomial in the DFT domain.
 * @param do_normalization The function normalizes the GLWE ciphertext if and only if do_normalization = 1.
 */
int const_mult_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw, const PolyUnivDFT* u,
                    int do_normalization);

//! GGSW IN DFT PART
// The same functions but the operation are done in the DFT domain. The parameters are the same but the ciphertext are in DFT
// space.

typedef struct ggsw_ciphertext_dft
{
	const GGSWCtParams* params;  // GGSW parameters
	MatBivDFT* mat;  // Represultent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in Zn[X]
} GGSWCiphertextDFT;

/**
 * @brief The number of coefficient in a bivariate GGSW ciphertext in the DFT domain.
 *
 * @param params_ggsw The GGSW parameters.
 * @return int64_t
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * Zn[X], due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t ggsw_coef_number_dft(const GGSWCtParams* params_ggsw);

/**
 * @brief Creates a bivGGSW in the DFT domain, filled with 0.
 *
 * @param params_ggsw The GGSW parameters.
 */
GGSWCiphertextDFT* new_ggsw_dft(const GGSWCtParams* params_ggsw);

/**
 * @brief Delete a bivGGSW in the DFT domain.
 *
 * @param ggsw The bivGGSW in the DFT domain.
 */
void delete_ggsw_dft(GGSWCiphertextDFT* ggsw_dft);

/**
 * @brief Return the pointer to biGLWE(DFT(-m * sk_j / 2^kappa*(i+1))) in the DFT domain.
 *
 * @param params_ggsw The GGSW parameters.
 * @param ggsw_mat The GGSW ciphertext's matrix.
 * @param j The j-th component of Sk.
 * @param i The i in the phase = -m * sk_j / (2^kappa_tilde)^(i+1).
 *
 * @return VecBivDFT*
 */
VecBivDFT* ggsw_Sj_Yti_dft(const GGSWCtParams* params_ggsw, MatBivDFT* ggsw_mat, int64_t j, int64_t i);

/**
 * @brief Normalize a GGSW ciphertext in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param result_dft The result normalized GGSW ciphertext in the DFT domain.
 * @param ggsw_dft The input GGSW ciphertext in the DFT domain.
 */
int normalize_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft);

/**
 * @brief Adds two GGSW ciphertexts in the DFT domain.
 *
 * @param result_dft The result GGSW ciphertext in the DFT domain.
 * @param ggsw_lhs_dft The left-hand side GGSW ciphertext in the DFT domain.
 * @param ggsw_rhs_dft The right-hand side GGSW ciphertext in the DFT domain.
 */
void add_ggsw_dft(GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_lhs_dft, const GGSWCiphertextDFT* ggsw_rhs_dft);

/**
 * @brief  Multiply a GGSW ciphertext by a constant in Zn[X]
 *
 * @param module Additionnal information for backend.
 * @param result_dft The result GGSW ciphertext.
 * @param ggsw_dft The GGSW ciphertext.
 * @param u_dft The Zn[X] polynomial in the DFT domain, with coefficient in [-2^(kappa-1), 2^(kappa-1)].
 * @param do_normalization The function normalizes the GLWE ciphertext if and only if do_normalization = 1.
 */
int const_mult_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft,
                        const PolyUnivDFT* u_dft, int do_normalization);

//! COMMON PART (begin)

/**
 * @brief Return the size of a bivGGSW ciphertext, in the DFT domain & out of the DFT domain.
 *
 * @param params_ggsw The GGSW parameters.
 * @return int64_t
 *
 * @note The size of a bivGGSW ciphertext is the same in and out of the DFT domain.
 */
uint64_t ggsw_size(const GGSWCtParams* params_ggsw);

/**
 * @brief The number of bytes needed to store a bivGGSW ciphertext.
 *
 * @param params_ggsw The GGSW parameters.
 * @return int64_t
 *
 * @note The number of bytes needed to store a bivGGSW ciphertext, is the same in and out of the DFT domain.
 */
uint64_t ggsw_bytes(const GGSWCtParams* params_ggsw);

#endif  // GGSW_CIPHERTEXT_H
