#ifndef GGSW_KEY_H
#define GGSW_KEY_H

#include "ggsw_ciphertext.h"

typedef struct ggsw_secret_key {
  uint64_t N;
  uint64_t k;
  PolyUniv** values;   // The key itself.
} GGSWSecretKey;

typedef struct ggsw_prep_secret_key {
  uint64_t N;
  uint64_t k;
  PolyUnivDFT** values;  // vec of size k, each element is prepared vec
  void* data;
} GGSWSecretKeyDFT;


//! GGSW KEY PART (begin)

/**
 * @brief Creates a new secret key values component.
 * 
 * @param N The dimension of the module Zn[X].
 * @param k The number of Zn[X] polynomials in a secret key.
 * @return PolyUniv** 
 */
PolyUniv** new_ggsw_secret_key_values(uint64_t N, uint64_t k);

/**
 * @brief Creates a new secret key values component, where each Zn[X] poly is uniformly drawn.
 * 
 * @param k The number of Zn[X] polynomials in a secret key.
 * @param N The dimension of the module Zn[X].
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 * @return PolyUniv** 
 */
PolyUniv** new_uniform_ggsw_secret_key_values(uint64_t N, uint64_t k, uint64_t  nb_bits);

/**
 * @brief Delete a secret key values component.
 * 
 * @param values The values of the secret key.
 * @param k The number of Zn[X] polynomials.
 */
void delete_ggsw_secret_key_values(PolyUniv** values, uint64_t k);

/**
 * @brief Creates a GGSW Secret key.
 * 
 * @param values The values of the secret key.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GGSWSecretKeyDFT* 
 */
GGSWSecretKey* new_ggsw_secret_key(PolyUniv** values, uint64_t N, uint64_t k);

/**
 * @brief Draws a secret key uniformly.
 * 
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 * @return GGSWSecretKeyDFT* 
 */
GGSWSecretKey* new_uniform_ggsw_secret_key(uint64_t N, uint64_t k, uint64_t  nb_bits);

/**
 * @brief Delete the secret key.
 * 
 * @param sk The secret key.
 */
void delete_ggsw_secret_key(GGSWSecretKey* sk);

/**
 * @brief Computes the secret key out of DFT space.
 * 
 * @param sk_dft The secret key in DFT space.
 * @return GGSWSecretKey* 
 */
GGSWSecretKey* transform_ggsw_secret_key_dft_to_not_dft(GGSWSecretKeyDFT* sk_dft);

/**
 * @brief Computes the values of the secret key out of DFT space.
 * 
 * @param values_dft The values of the secret key in DFT space.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return PolyUniv** 
 */
PolyUniv** transform_ggsw_secret_key_values_dft_to_not_dft(PolyUnivDFT** values_dft, uint64_t N, uint64_t k);


//! GGSW IN DFT SPACE PART (begin)

/**
 * @brief Creates a new secret key values component in DFT space.
 * 
 * @param N The dimension of the module Zn[X].
 * @param k The number of Zn[X] polynomials in a secret key.
 * @return PolyUnivDFT** 
 */
PolyUnivDFT** new_ggsw_secret_key_values_dft(uint64_t N, uint64_t k);

/**
 * @brief Creates a new secret key values component in DFT space, where each Zn[X] poly is uniformly drawn.
 * 
 * @param k The number of Zn[X] polynomials in a secret key.
 * @param N The dimension of the module Zn[X].
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 * @return PolyUnivDFT** 
 */
PolyUnivDFT** new_uniform_ggsw_secret_key_values_dft(uint64_t N, uint64_t k, uint64_t  nb_bits);

/**
 * @brief Delete the values of a secret key in DFT space.
 * 
 * @param values The values of the secret key in DFT space.
 */
void delete_ggsw_secret_key_values_dft(PolyUnivDFT** values, uint64_t k);

/**
 * @brief Creates a GGSW Secret key in DFT space.
 * 
 * @param values The values of the secret key in DFT space.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GGSWSecretKeyDFT* 
 */
GGSWSecretKeyDFT* new_ggsw_secret_key_dft(PolyUnivDFT** values, uint64_t N, uint64_t k);

/**
 * @brief Draws a secret key uniformly in DFT space.
 * 
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GGSWSecretKeyDFT* 
 */
GGSWSecretKeyDFT* new_uniform_ggsw_secret_key_dft(uint64_t N, uint64_t k, uint64_t nb_bits);

/**
 * @brief Delete the secret key that is in DFT space.
 * 
 * @param sk_dft The secret key in DFT space.
 */
void delete_ggsw_secret_key_dft(GGSWSecretKeyDFT* sk_dft);

/**
 * @brief Computes the secret key in DFT space.
 * 
 * @param sk The secret key out of DFT space.
 * @return GGSWSecretKey* 
 */
GGSWSecretKeyDFT* transform_ggsw_secret_key_not_dft_to_dft(GGSWSecretKey* sk);

/**
 * @brief Computes the values of the secret key in DFT space.
 * 
 * @param values The values of the secret key out of DFT space.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return PolyUnivDFT** 
 */
PolyUnivDFT** transform_ggsw_secret_key_values_not_dft_to_dft(PolyUniv** values, uint64_t N, uint64_t k);

typedef struct ggsw_public_key {
  uint32_t size;       // number of ciphertexts in public key
  GGSWCiphertext* pk;  // Public key is multiple encryptions of 0
} GGSWPublicKey;

#endif  // GGSW_KEY_H
