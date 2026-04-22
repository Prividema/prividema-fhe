#ifndef bivGGSW_CT_PARAMS_H
#define bivGGSW_CT_PARAMS_H

#include "glwe_params.h"

/**
 * @struct GGSWParams
 *
 * @brief GGSW parameters.
 */
typedef struct ggsw_ct_params
{
	/// GLWE Parameters.
	const GLWEParams* params_glwe;
	/**
	 * @brief \f$ \tilde{k} \f$ parameter
	 *
	 * By default \f$ k = \tilde{k} \f$, which for RGSW is = 1
	 */
	uint64_t k_tilde;
	/**
	 * @brief \kappa_tilde, the exponent that will be used to multiply the message by (negative) powers of 2 thereof.
	 *
	 * Imporntantly, \K and \kappa_tilde are different parameters of a GGSW Ciphertext, one encodes the radix size and the other
	 * the powers used to multiply the message repeatedly.
	 */
	uint64_t kappa_tilde;
	/**
	 * @brief Number of total "limbs" (bivariate GLWE ciphertexts) in the GGSWParams.
	 *
	 * This is, therefore,
	 * \f$ \tilde{l}_a \cdot k + \tilde{l}_b \f$
	 *
	 *
	 * As in the GLWE case, only \f$ \tilde{l}_a = \tilde{l}_b \f$ and \f$ \tilde{l}_a = \tilde{l}_b + 1 \f$
	 * are supported
	 */
	uint64_t ciphertext_nb_limbs_tilde;
} GGSWParams;

/**
 * @brief Creates a set of GGSW Parameters.
 *
 * @param params         The GLWE parameters.
 * @param k_tilde        \f$ \tilde{k} \f$ the number of polynomials in the secret key.
 * @param kappa_tilde    \kappa_tilde
 * @param nb_limbs_tilde \f$ \tilde{l}_a \cdot k + \tilde{l}_b \f$
 *
 * @return A Pointer to the set of GGSW Parameters, or NULL in case of failure
 */
GGSWParams* new_ggsw_params(const GLWEParams* params, uint64_t k_tilde, uint64_t kappa_tilde, uint64_t nb_limbs_tilde);

/**
 * @brief Deletes a set of GGSW parameters.
 *
 * @param params A Pointer to the GGSW parameters to delete.
 */
void delete_ggsw_params(GGSWParams* params);

/**
 * @brief Gets the number of rows in a GGSW.
 *  In other words, the number of GLWEs that form a GGSW,
 *  \f$ \tilde{l}_a \cdot k + \tilde{l}_b \f$
 *
 * @param params The GGSW parameters.
 *
 * @return The number of rows in a GGSW.
 */
uint64_t ggsw_num_rows(const GGSWParams* params);

/**
 * @brief Gets the parameter \f$ \tilde{l}_a \f$
 *
 * @param params The GGSW parameters
 *
 * @return \f$ \tilde{l}_a \f$
 */
uint64_t ggsw_params_l_tilde_a(const GGSWParams* params);

/**
 * @brief Gets the parameter \f$ \tilde{l}_b \f$
 *
 * @param params The GGSW parameters
 *
 * @return \f$ \tilde{l}_b \f$
 */
uint64_t ggsw_params_l_tilde_b(const GGSWParams* params);

/**
 * @brief Computes the number of coefficients in a Bivariate GGSW ciphertext.
 *
 * @param params_ggsw The GGSW parameters.
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
 * @remark The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * \ZnX, due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t ggsw_coef_number_dft(const GGSWParams* params_ggsw);

/**
 * @brief Gets the number of bytes needed to store a GGSW ciphertext.
 *
 * @param params_ggsw The GGSW parameters.
 *
 * @return The number of bytes needed to store a GGSW ciphertext.
 *
 * @remark The result is the same for a coefficient-domain or DFT domain representation
 */
uint64_t ggsw_bytes(const GGSWParams* params_ggsw);

/**
 * @brief Gets the total number of polynomials that make up a GGSW
 *
 * Named because each polynomial is a "limb" of a GLWE, of which there are
 * num_rows of in a GGSW.
 *
 * @param params_ggsw The GGSW parameters.
 *
 * @return The size (in terms of polynomials) of a GGSW ciphertext..
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
	uint64_t kappa_tilde;  ///< \kappa_tilde as in the papers
	uint64_t l_tilde;      ///< The precision, \f$ \tilde{l} \f$ as in the papers
} GLWEGadgetParams;

/**
 * @brief Creates a set of GGSW Parameters.
 *
 * @param params        The GLWE parameters.
 * @param kappa_tilde   The 2-exponent of the base (\f$ 2^{-\kappatilde} \f$)
 * @param l_tilde       \f$ \tilde{l} \f$
 *
 * @return A Pointer to the set of GGSW Parameters, or NULL if it failed to create it.
 */
GLWEGadgetParams* new_glwegadget_params(const GLWEParams* params, uint64_t kappa_tilde, uint64_t l_tilde);

/**
 * @brief Deletes a glwegadget's parameters
 *
 * @param params The GLWEGadgetParams to delete/deallocate
 */
void* delete_glwegadget_params(GLWEGadgetParams* params);

#endif  // bivGGSW_CT_PARAMS_H
