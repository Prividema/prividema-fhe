#ifndef GLWEGADGET_KEY_H
#define GLWEGADGET_KEY_H

#include <stdint.h>

#include "glwegadget_ciphertext.h"

/**
 * @file glwegadget_key.h
 *
 * @brief Contains data structures related GLWEGadget specific keys, mainly key-switching-keys.
 *
 * So far only KSK for automprphisms are here
 */

typedef struct glwe_automorphism_keyswitching_key_t
{
	int automorphism_p;                ///< -2 represents uninitialized value, 0 if the KSK is not for an automorphism
	GLWEGadgetParams* params;          ///< The GLWEGadgetParams associated with the gadgets used for key switching
	GLWEGadgetCiphertextPrep** enc_s;  ///< A vector with pointers to the k prepared ciphertetxts of f(sk_i)
} GLWEAutomorphismKSK;

/**
 * @brief Creates a new Key-switching-key (KSK)
 *
 * @param params The GLWEGadgetParams associated with the KSK
 *
 * @return NULL in case of failure, a new KSK otherwise
 *
 */
GLWEAutomorphismKSK* new_automorphism_ksk(GLWEGadgetParams* params);

/**
 * @brief Deletes a KSK
 *
 * @param automorphism_ksk The KSK to delete
 */
void delete_automorphism_ksk(GLWEAutomorphismKSK* automorphism_ksk);

/**
 * Collection of KSKs
 *
 */
typedef struct glwe_automorphism_keyswitching_key_collection_t
{
	uint64_t size;               ///< Number of positions in the keys vector
	GLWEAutomorphismKSK** keys;  ///< Pointers to KSKs
} GLWEAutomorphismKSKCollection;

/**
 * @brief Creates a KSK collection
 *
 * @param size The number of KSKs that the collection may contain
 *
 * @return A new empty KSK collection, with no KSK present nor allocated. NULL in case of failure
 */
GLWEAutomorphismKSKCollection* new_automorphism_ksk_collection(uint64_t size);

/**
 * @brief Puts a KSK into the specified position in the collection. Returns previous occupant.
 *
 * @param collection The KSK collection to put the KSK in
 * @param key The KSK to put (a pointer of) in the collection
 * @param pos The position in the collection where the KSK will be
 *
 * @return (a pointer to) The KSK that previously occupied position pos, NULL if unoccupied
 */
GLWEAutomorphismKSK* glwegadget_ksk_collection_put_key(GLWEAutomorphismKSKCollection* collection,
                                                       GLWEAutomorphismKSK* key, uint64_t pos);

/**
 * @brief Retrieves the KSK in position pos from a collection
 *
 * @param collection The KSK collection to fetch the KSK from
 * @param pos The position in the collection to fetch
 *
 * @return (a pointer to) The KSK  occupies position pos, NULL if unoccupied
 */
GLWEAutomorphismKSK* glwegadget_ksk_collection_get_key(const GLWEAutomorphismKSKCollection* collection, uint64_t pos);

/**
 * @brief Deletes a KSK collection, optionally deallocating the KSKs in it
 *
 * @param automorphism_ksk The collection to delete
 * @param deallocate_ksks Whether to delete the KSKs pointed by it (1) or not (0)
 *
 * @return The number of KSKs that have been deleted while deleting the collection
 */
uint64_t delete_automorphism_ksk_collection(GLWEAutomorphismKSKCollection* automorphism_ksk, int deallocate_ksks);

#endif
