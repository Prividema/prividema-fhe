#ifndef GGSW_CIPHERTEXT_H
#define GGSW_CIPHERTEXT_H

/**
 * @file ggsw_ciphertext.h
 * @brief bivGGSW ciphertext
 *
 * In this header file, we define the structure represultenting bivariate bivGGSW ciphertext
 * and bivariate bivGGSW in the DFT domain. Moreover, we define alias for double and int64_t.
 * Thus, we properly manipulate mathematical structure (polynomial, vector or matrix).
 * The types below are provided solely for clarity, nothing is changed in spqlios.
 */

#include <stdint.h>
#include <stdlib.h>

#include "ggsw_ct_params.h"
#include "maths_structures.h"
#include "partialggsw_ciphertext.h"

//! bivGGSW PART (begin)

typedef struct ggsw_ciphertext
{
	const GGSWCtParams* params;  // bivGGSW parameters
	MatBiv* mat;  // Represultent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in Zn[X]
} GGSWCiphertext;

/**
 * @brief Return the number of coefficients in a bivariate bivGGSW ciphertext.
 *
 * @param params_ggsw The bivGGSW parameters.
 * @return int64_t
 */
uint64_t ggsw_coef_number(const GGSWCtParams* params_ggsw);

/**
 * @brief Creates a bivGGSW ciphertext, it is filled with 0.
 *
 * @param params_ggsw The bivGGSW parameters.
 *
 * @retval - `NULL` if malloc failed inside the function.
 * @retval - Allocated bivGGSW ciphertext, filled with 0.
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
 * @brief Return the pointer to : biGLWE(DFT(-m * sk_j / 2^kappa_tilde^i)) if j < k.
 * 								  biGLWE(DFT(m / 2^kappa_tilde^i))         if j = k.
 *
 * @param params_ggsw The bivGGSW parameters.
 * @param ggsw_mat The bivGGSW ciphertext's matrix.
 * @param j The j-th component of the secret key.
 * @param i The i in -m * sk_j / (2^kappa_tilde)^i.
 *
 * @return VecBiv*
 */
VecBiv* ggsw_Sj_Yti(const GGSWCtParams* params_ggsw, MatBiv* ggsw_mat, int64_t j, int64_t i);

/**
 * @brief Normalize a bivGGSW ciphertext.
 *
 * @param module Additionnal information for backend.
 * @param result The result normalized bivGGSW ciphertext.
 * @param ggsw The input bivGGSW ciphertext.
 * 
 * @retval - `-1` if malloc failed inside the function vec_znx_normalize_base2k_p.
 * @retval - `0`.
 */
int normalize_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw);

/**
 * @brief Adds two bivGGSW ciphertexts with same params and put result in result
 *
 * @param result The result bivGGSW ciphertext.
 * @param ggsw_lhs The left-hand side bivGGSW ciphertext.
 * @param ggsw_rhs The right-hand side bivGGSW ciphertext.
 */
void add_ggsw(GGSWCiphertext* result, const GGSWCiphertext* ggsw_lhs, const GGSWCiphertext* ggsw_rhs);

/**
 * @brief Multiply a bivGGSW ciphertext by a constant in Zn[X]
 *
 * @param module Additionnal information for backend.
 * @param result The result bivGGSW ciphertext.
 * @param ggsw The bivGGSW ciphertext.
 * @param u_dft The Zn[X] polynomial in the DFT domain.
 * @param do_normalization The function normalizes the GLWE ciphertext if and only if do_normalization = 1.
 */
int const_mult_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw, const PolyUnivDFT* u_dft,
                    int do_normalization);

//! bivGGSW IN DFT PART
// The same functions but the operation are done in the DFT domain. The parameters are the same but the ciphertext are in the DFT 
// domain.

typedef struct ggsw_ciphertext_dft
{
	const GGSWCtParams* params;  // bivGGSW parameters
	MatBivDFT* mat;  // Represent a matrix of size n_limbs_tilde x n_limbs with coefficients that are in Cn/2[X]
} GGSWCiphertextDFT;

/**
 * @brief The number of coefficient in a bivariate bivGGSW ciphertext in the DFT domain.
 *
 * @param params_ggsw The bivGGSW parameters.
 * @return int64_t
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * Zn[X], due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t ggsw_coef_number_dft(const GGSWCtParams* params_ggsw);

/**
 * @brief Creates a bivGGSW in the DFT domain, filled with 0.
 *
 * @param params_ggsw The bivGGSW parameters.
 * 
 * @retval - `NULL` if malloc failed inside the function.
 * @retval - Allocated bivGGSW ciphertext in the DFT domain, filled with 0.
 */
GGSWCiphertextDFT* new_ggsw_dft(const GGSWCtParams* params_ggsw);

/**
 * @brief Delete a bivGGSW in the DFT domain.
 *
 * @param ggsw The bivGGSW in the DFT domain.
 */
void delete_ggsw_dft(GGSWCiphertextDFT* ggsw_dft);

/**
 * @brief Return the pointer to : biGLWE(DFT(-m * sk_j / 2^kappa_tilde^i)) if j < k in the DFT domain.
 * 								  biGLWE(DFT(m / 2^kappa_tilde^i))         if j = k in the DFT domain.
 *
 * @param params_ggsw The bivGGSW parameters.
 * @param ggsw_mat_dft The bivGGSW ciphertext's matrix in the DFT domain.
 * @param j The j-th component of Sk.
 * @param i The i in the phase = -m * sk_j / (2^kappa_tilde)^i.
 *
 * @return VecBivDFT*
 */
VecBivDFT* ggsw_Sj_Yti_dft(const GGSWCtParams* params_ggsw, MatBivDFT* ggsw_mat_dft, int64_t j, int64_t i);

/**
 * @brief Normalize a bivGGSW ciphertext in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param result_dft The result normalized bivGGSW ciphertext in the DFT domain.
 * @param ggsw_dft The input bivGGSW ciphertext in the DFT domain.
 * 
 * @retval `-1` if malloc failed inside vec_znx_normalize_base2k_p.
 * @retval `0`.
 */
int normalize_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft);

/**
 * @brief Adds two bivGGSW ciphertexts in the DFT domain.
 *
 * @param result_dft The result bivGGSW ciphertext in the DFT domain.
 * @param ggsw_lhs_dft The left-hand side bivGGSW ciphertext in the DFT domain.
 * @param ggsw_rhs_dft The right-hand side bivGGSW ciphertext in the DFT domain.
 */
void add_ggsw_dft(GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_lhs_dft, const GGSWCiphertextDFT* ggsw_rhs_dft);

/**
 * @brief Multiply a bivGGSW ciphertext by a constant in Zn[X].
 *
 * @param module Additionnal information for backend.
 * @param result_dft The result bivGGSW ciphertext in the DFT domain.
 * @param ggsw_dft The bivGGSW ciphertext in the DFT domain.
 * @param u_dft The Zn[X] polynomial in the DFT domain, with coefficient in [-2^(kappa-1), 2^(kappa-1)].
 * @param do_normalization The function normalizes the GLWE ciphertext if and only if do_normalization = 1.
 * 
 * @retval `-1` if malloc failed inside the function.
 * @retval `0`.
 */
int const_mult_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft,
                        const PolyUnivDFT* u_dft, int do_normalization);

//! COMMON PART (begin)

/**
 * @brief Return the size of a bivGGSW ciphertext, in the DFT domain & out of the DFT domain.
 *
 * @param params_ggsw The bivGGSW parameters.
 * @return int64_t
 *
 * @note The size of a bivGGSW ciphertext is the same in and out of the DFT domain.
 */
uint64_t ggsw_size(const GGSWCtParams* params_ggsw);

/**
 * @brief The number of bytes needed to store a bivGGSW ciphertext.
 *
 * @param params_ggsw The bivGGSW parameters.
 * @return int64_t
 *
 * @note The number of bytes needed to store a bivGGSW ciphertext, is the same in and out of the DFT domain.
 */
uint64_t ggsw_bytes(const GGSWCtParams* params_ggsw);

#endif  // GGSW_CIPHERTEXT_H
