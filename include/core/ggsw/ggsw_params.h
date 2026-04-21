#ifndef bivGGSW_CT_PARAMS_H
#define bivGGSW_CT_PARAMS_H

/**
 * @file ggsw_ct_params.h
 *
 * In this header, we define the structure representing bivGGSW parameters.
 */

#include "glwe_params.h"

/**
 * @struct GGSWCtParams
 *
 * @brief GGSW parameters.
 */
typedef struct ggsw_ct_params
{
	/// GLWE Parameters.
	const GLWEParams* params_glwe;
	/**
	 * By default k = k_tilde, which for RGSW is = 1
	 */
	uint64_t k_tilde;
	/**
	 * kappa tilde, the exponent that will be used to multiply the message by (negative) powers of 2 thereof.
	 *
	 * Imporntantly, kappa and kappa_tilde are different parameters of a GGSW Ciphertext, one encodes the radix size and the other
	 * the powers used to multiply the message repeatedly.
	 */
	uint64_t kappa_tilde;
	/**
	 *Number of total limbs in the GGSWParams. This is, therefore, l_tilde_a * k_tilde + l_tilde_b
	 */
	uint64_t ciphertext_nb_limbs_tilde;
} GGSWParams;

/**
 * @brief Creates a set of GGSW Parameters.
 *
 * @param params_glwe    The GLWE parameters.
 * @param k_tilde        The number of polynomials in the secret key.
 * @param kappa_tilde    Kappa tilde
 * @param nb_limbs_tilde l_tilde_a * k_tilde + l_tilde_b
 *
 * @return A Pointer to the set of GGSW Parameters.
 */
GGSWParams* new_ggsw_params(const GLWEParams* params, uint64_t k_tilde, uint64_t kappa_tilde, uint64_t nb_limbs_tilde);

/**
 * @brief Deletes a GGSW parameters.
 *
 * @param params A Pointer to the GGSW parameters.
 */
void delete_ggsw_params(GGSWParams* params);

/**
 * @brief Gets the number of rows in a GGSW.
 *  In other words, the number of GLWEs that form a GGSW,
 *  or l_tilde_a * k_tilde + l_tilde_b
 *
 * @param params The GGSW parameters.
 *
 * @return The number of rows in a GGSW.
 */
uint64_t ggsw_num_rows(const GGSWParams* params);

/**
 * @brief Gets the parameter l_tilde_a
 *
 * @param params The GGSW params
 *
 * @return l_tilde_a
 */
uint64_t ggsw_params_l_tilde_a(const GGSWParams* params);

/**
 * @brief Gets the parameter l_tilde_b
 *
 * @param params The GGSW params
 *
 * @return l_tilde_b
 */
uint64_t ggsw_params_l_tilde_b(const GGSWParams* params);

/**
 * @brief Computes the number of coefficients in a Bivariate GGSW ciphertext.
 *
 * @param params_ggsw A Pointer to the GGSW parameters.
 *
 * @return The number of coefficients in a GGSW ciphertext.
 */
uint64_t ggsw_coef_number(const GGSWParams* params_ggsw);

/**
 * @brief Computes the number of coefficients in a Bivariate GGSW ciphertext in the DFT space.
 *
 * @param params_ggsw A Pointer to the GGSW parameters.
 *
 * @return The number of coefficients in a GGSW ciphertext.
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * \f$\mathbb{Z}_n[X]\f$, due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t ggsw_coef_number_dft(const GGSWParams* params_ggsw);

/**
 * @brief Gets the number of bytes needed to store a GGSW ciphertext.
 *
 * @param params_ggsw A Pointer to the GGSW parameters.
 *
 * @return The number of bytes needed to store a GGSW ciphertext.
 *
 * @note This function works in both DFT and iDFT domains.
 */
uint64_t ggsw_bytes(const GGSWParams* params_ggsw);

/**
 * @brief Gets the size of any type of GGSW ciphertext.
 *
 * @param params_ggsw A Pointer to the GGSW parameters.
 *
 * @return The size of any type of GGSW ciphertext.
 *
 * @note The size is the same in DFT and iDFT domains.
 */
uint64_t ggsw_total_n_glwe_limbs(const GGSWParams* params_ggsw);

/**
 * @struct GLWEGadgetParams
 *
 * @brief GLWEGadget Parameters.
 */
typedef struct glwegadget_params
{
	/// GLWE parameters.
	const GLWEParams* params_glwe;
	uint64_t kappa_tilde;  ///< Kappa tilde as in the papers
	uint64_t l_tilde;      ///< The precision, l_tilde as in the papers
} GLWEGadgetParams;

/**
 * @brief Creates a set of GGSW Parameters.
 *
 * @param params_glwe 	The GLWE parameters.
 * @param kappa_tilde 	The 2-exponent of the base (2^kappa_tilde).
 * @param l_tilde
 *
 * @return A Pointer to the set of GGSW Parameters.
 */
GLWEGadgetParams* new_glwegadget_params(const GLWEParams* params, uint64_t kappa_tilde, uint64_t l_tilde);

/**
 * @brief Deletes a glwegadget's parameters
 *
 * @param params The GLWEGadgetParams to delete/deallocate
 */
void* delete_glwegadget_params(GLWEGadgetParams* params);

#endif  // bivGGSW_CT_PARAMS_H
