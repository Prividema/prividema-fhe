#ifndef PARTIALGGSW_CIPHERTEXT_H
#define PARTIALGGSW_CIPHERTEXT_H

#include <stdint.h>

#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "glwe_key.h"

/**
 *
 * @file glwegadget_ciphertext.h
 *
 * @brief Definition of a GLWEGadget ciphertext, creation, destruction, encryption and transformations thereof
 *
 * Decrytion is not implemented nor is it planned:
 * GLWEGadgets are used as an input format for half-external products, not as the main format for encrypted data.
 * They contain multiple encryptions of the same message multiplied by different constants.
 * Thus, decrypting them does not seem to be a useful operation, nor a clearly-defined one:
 * should we output all the products of the message by powers of \f$ 2^{\tilde{K}} \f$ ?
 * Should we check that the decryption correctly gives us such a sequence of plaintexts?
 *
 * Given all of that, we decided not to provide a decryption function.
 * If needed (as in the tests, which one can take as an example) the GLWE decryption function(s) can be used
 * together with the glwegadget_extract_bivglwe function, but, again, outside of debugging it is
 * higly unlikely one might need to to that.
 */

/**
 * @brief A GLWEGadget ciphertext
 */
typedef struct glwegadget_ciphertext
{
	/// GLWEGadgetParams (parameters) of this ciphertext
	const GLWEGadgetParams* params;

	/// @brief internal GLWEGadget data
	/// GLWEGadget data stored as a flattened row-major
	/// matrix of size l_tilde * n_limbs, each row the contents of a GLWECiphertext
	///
	/// See \ref glwegadget_encoding "GLWEGadget's memory layout expanation" for more details
	MatBiv* mat;
} GLWEGadgetCiphertext;

/**
 * @brief creates a new GLWEGadgetCiphertext
 *
 * @param params The parameters
 *
 */
GLWEGadgetCiphertext* new_glwegadget(const GLWEGadgetParams* params);

/**
 * @brief Deletes a GLWEGadget ciphertext
 *
 * @param glwegadget_ct The ciphertext to delete
 *
 */
void delete_glwegadget(GLWEGadgetCiphertext* glwegadget_ct);

/**
 * @brief Debugging function to inspect GLWEGadgets
 *
 * Intended to inspect GLWEGadget contents for debugging.
 * It takes as input a secret key to decrypt. For debugging, a global variable could
 * be used to store the reference.
 *
 * Prints the n first (lowest degree) coefficients of the decryption.
 * Only able to output the 64 MSBs of the encrypted message/phase
 *
 * @param module The backend module
 * @param glwe_gad The GLWEGadget to print
 * @param sk_prep The prepared secret key for the GLWEGadget
 * @param n How many coefficients to print
 *
 */
int print_coefs_gad(const PvdaBackend* module, const GLWEGadgetCiphertext* glwe_gad,
                    const GLWESecretKeyPrepared* sk_prep, int n);

/**
 * @brief A GLWEGadget ciphertext that has been preprocessed for use in a half-external product
 */
typedef struct glwegadget_ciphertext_prepared
{
	/// GLWEGadgetParams (parameters) of this ciphertext
	const GLWEGadgetParams* params;
	MatBivDFT*
	    mat;  ///< Prepared GLWEGadgetCiphertext matrix. Its layout is an internal implementation detail of the underlying compute engine
} GLWEGadgetCiphertextPrep;

/**
 * @brief creates a new prepared GLWEGadget
 *
 * @param params The parameters
 *
 */
GLWEGadgetCiphertextPrep* new_glwegadget_prep(const GLWEGadgetParams* params);

/**
 * @brief Deletes a prepared GLWEGadget
 *
 * @param glwegadget_prep_ct The ciphertext to delete
 *
 */
void delete_glwegadget_prep(GLWEGadgetCiphertextPrep* glwegadget_prep_ct);

/**
 * @brief Encrypts a univariate integer polynomial into a GLWEGadget
 * (also known as halfGGSW in the k=1 case)
 *
 * @param module The underlying compute module
 * @param result The resulting GLWEGadget ciphertext
 * @param sk_prep A prepared secret key
 * @param m_univ The univariate plaintext
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwegadget_secret_encrypt(const PvdaBackend* module, GLWEGadgetCiphertext* result,
                              const GLWESecretKeyPrepared* sk_prep, const PolyUniv* m_univ);

/**
 * @brief Encrypts a univariate integer polynomial into a GLWEGadget
 * (also known as halfGGSW in the k=1 case)
 *
 *  TODO: not yet implemented
 *
 * @param module The underlying compute module
 * @param result The resulting GLWEGadget ciphertext
 * @param pk     A public key
 * @param m_univ The univariate plaintext
 *
 */
void glwegadget_public_encrypt(const PvdaBackend* module, GLWEGadgetCiphertext* result, const GLWEPublicKey* pk,
                               const PolyUniv* m_univ);
/**
 *
 * Encrypts a univariate polynomial with at most d non-zero coefficients
 * into a "packed GLWEGadget" GLWE ciphertext.
 *
 * See https://github.com/Prividema/prividema-fhe/pull/64 for an explanation on "packed" GLWEGadgets
 *
 * TODO: move documentation to the docs folder instead of the PR
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwegadget_packed_secret_encrypt(const PvdaBackend* module, GLWECiphertext* result,
                                     const GLWEGadgetParams* params_glwegad, const GLWESecretKeyPrepared* sk_prep,
                                     const PolyUniv* m_univ, uint64_t d);

/**
 * @brief Retreive a bivglwe from a GLWEGadget (halfGGSW for k = 1)
 *
 * @param glwegadget_ct     The GLWEGadget ciphertext
 * @param i                 The index of the GLWE to retrievie, from 1 to l_tilde included
 *
 */
VecBiv* glwegadget_extract_bivglwe(GLWEGadgetCiphertext* glwegadget_ct, uint64_t i);

/**
 * @brief Prepares a GLWEGadget for use in a half-external product
 *
 * @param module The underlying backend module
 * @param glwegadget_prep_ct The resulting prepared GLWEGadget
 * @param glwegad_ct The input GLWEGadget
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 *
 */
int glwegadget_prepare(const PvdaBackend* module, GLWEGadgetCiphertextPrep* glwegadget_prep_ct,
                       const GLWEGadgetCiphertext* glwegad_ct);

#endif  // PARTIALGGSW_CIPHERTEXT_H
