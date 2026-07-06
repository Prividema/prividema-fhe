#ifndef GLWEGADGET_KEY_H
#define GLWEGADGET_KEY_H

#include <stdint.h>

#include "glwegadget_ciphertext.h"

/**
 * @file glwegadget_key.h
 *
 * @brief Contains data structures related GLWEGadget specific keys, mainly key-switching-keys
 * (and the automorphism keys).
 *
 */

typedef struct glwe_ksk_t
{
	GLWEGadgetParams* params;          ///< The GLWEGadgetParams associated with the gadgets used for key switching
	GLWEGadgetCiphertextPrep** enc_s;  ///< A vector with pointers to the k prepared ciphertetxts of f(sk_i)
} GLWEKSK;

/**
 * Create a GLWE Key-switching key (KSK)
 *
 * @param params The GLWEGadgetParams associated with the KSK
 *
 * @return NULL in case of failure, a new KSK otherwise
 *
 */
GLWEKSK* new_glwe_ksk(GLWEGadgetParams* params);
/**
 * @brief Deletes a GLWE KSK
 *
 * @param glwe_ksk The GLWE KSK to delete
 */
void delete_glwe_ksk(GLWEKSK* glwe_ksk);

typedef struct glwe_automorphism_keyswitching_key_t
{
	int automorphism_p;  ///< -2 represents uninitialized value, 0 if the KSK is not for an automorphism
	GLWEKSK* glwe_ksk;
} GLWEAutomorphismKey;

/**
 * @brief Creates a new automorphism key
 *
 * @param params The GLWEGadgetParams associated with the automorphism key
 *
 * @return NULL in case of failure, a new automorphism key otherwise
 *
 */
GLWEAutomorphismKey* new_automorphism_key(GLWEGadgetParams* params);

/**
 * @brief Deletes an automorphism key
 *
 * @param automorphism_key The automorphism key to delete
 */
void delete_automorphism_key(GLWEAutomorphismKey* automorphism_key);

/**
 * Collection of automorphism keys
 *
 */
typedef struct glwe_automorphism_key_collection_t
{
	uint64_t size;               ///< Number of positions in the keys vector
	GLWEAutomorphismKey** keys;  ///< Pointers to automorphism keys
} GLWEAutomorphismKeyCollection;

/**
 * @brief Creates an automorphism key collection
 *
 * @param size The number of automorphism keys that the collection may contain
 *
 * @return A new empty collection, with no automorphism keys present nor allocated. NULL in case of failure
 */
GLWEAutomorphismKeyCollection* new_automorphism_key_collection(uint64_t size);

/**
 * @brief Puts an automorphism key into the specified position in the collection. Returns previous occupant.
 *
 * @param collection The automorphism key collection to put the automorphism key in
 * @param key The automorphism key to put (a pointer of) in the collection
 * @param pos The position in the collection where the automorphism key will be
 *
 * @return (a pointer to) The automorphism key that previously occupied position pos, NULL if unoccupied
 */
GLWEAutomorphismKey* glwegadget_key_collection_put_key(GLWEAutomorphismKeyCollection* collection,
                                                       GLWEAutomorphismKey* key, uint64_t pos);

/**
 * @brief Retrieves the automorphism key in position pos from a collection
 *
 * @param collection The collection to fetch the automorphism key from
 * @param pos The position in the collection to fetch
 *
 * @return (a pointer to) The automorphism occupying position pos, NULL if unoccupied
 */
GLWEAutomorphismKey* glwegadget_key_collection_get_key(const GLWEAutomorphismKeyCollection* collection, uint64_t pos);

/**
 * @brief Deletes an automorphism key collection, optionally deallocating the keys in it
 *
 * @param automorphism_keys The collection to delete
 * @param deallocate_keys Whether to delete the keys pointed by it (1) or not (0)
 *
 * @return The number of keys that have been deleted while deleting the collection
 */
uint64_t delete_automorphism_key_collection(GLWEAutomorphismKeyCollection* automorphism_keys, int deallocate_keys);

#endif
