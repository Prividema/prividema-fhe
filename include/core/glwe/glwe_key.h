#ifndef GLWE_KEY_H
#define GLWE_KEY_H

#include "glwe_ciphertext.h"

// TODO : To be Defined (spqlios).
/**
 * @brief obtain a info for:
 *  - the dimension N (or the complex dimension m=N/2)
 *  - any moduleuted fft or ntt items
 *  - the hardware (avx, arm64, x86, ...)
 */
typedef struct core {

} Core;

//! GLWE SECRET KEY STRUCTURES
typedef struct glwe_secret_key {
  uint64_t N;
  uint64_t k;
  PolyUniv** values;
  void* data;
} GLWESecretKey;

typedef struct glwe_prep_secret_key {
  uint64_t N;
  uint64_t k;
  PolyUnivDFT** values;  // vec of size k, each element is prepared vec
  void* data;
} GLWEPreparedSK;


//! GLWE PART (begin)

/**
 * @brief Creates a new secret key values component.
 * 
 * @param k The number of Zn[X] polynomials in a secret key.
 * @param N The dimension of the module Zn[X].
 * @return PolyUniv** 
 */
PolyUniv** new_secret_key_values(int64_t N, int64_t k);

/**
 * @brief Delete a secret key values component.
 * 
 * @param values The values of the secret key.
 * @param k The number of Zn[X] polynomials.
 */
void delete_secret_key_values(PolyUniv** values, int64_t k);


/**
 * @brief Creates a GLWE Secret key.
 * 
 * @param values The values of the secret key.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GLWEPreparedSK* 
 */
GLWESecretKey* new_glwe_secret_key(PolyUniv** values, int64_t N, int64_t k);

/**
 * @brief Delete the secret key.
 * 
 * @param sk_dft The secret key.
 */
void delete_glwe_secret_key(GLWESecretKey* sk);

/**
 * @brief Computes the secret key out of DFT space.
 * 
 * @param sk_dft The secret key in DFT space.
 * @return GLWESecretKey* 
 */
GLWESecretKey* transform_glwe_secret_key_dft_to_not_dft(GLWEPreparedSK* sk_dft);

/**
 * @brief Computes the values of the secret key out of DFT space.
 * 
 * @param sk_dft The secret key in DFT space.
 * @return PolyUniv** 
 */
PolyUniv** transform_secret_key_values_dft_to_not_dft(GLWEPreparedSK* sk_dft);


//! GLWE IN DFT SPACE PART (begin)

/**
 * @brief Creates a new secret key values component in DFT space.
 * 
 * @param k The number of Zn[X] polynomials in a secret key.
 * @param N The dimension of the module Zn[X].
 * @return PolyUnivDFT** 
 */
PolyUnivDFT** new_secret_key_values_dft(int64_t N, int64_t k);

/**
 * @brief Delete the values of a secret key in DFT space.
 * 
 * @param values The values of the secret key in DFT space.
 */
void delete_secret_key_values_dft(PolyUnivDFT** values, int64_t k);

/**
 * @brief Creates a GLWE Secret key in DFT space.
 * 
 * @param values The values of the secret key in DFT space.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GLWEPreparedSK* 
 */
GLWEPreparedSK* new_glwe_secret_key_dft(PolyUnivDFT** values, int64_t N, int64_t k);

/**
 * @brief Delete the secret key that is in DFT space.
 * 
 * @param sk_dft The secret key in DFT space.
 */
void delete_glwe_secret_key_dft(GLWEPreparedSK* sk_dft);

/**
 * @brief Computes the secret key in DFT space.
 * 
 * @param sk_dft The secret key out of DFT space.
 * @return GLWESecretKey* 
 */
GLWEPreparedSK* transform_glwe_secret_key_not_dft_to_dft(GLWESecretKey* sk);

/**
 * @brief Draws a secret key uniformly in DFT space.
 * 
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GLWEPreparedSK* 
 */
GLWEPreparedSK* new_uniform_glwe_secret_key_dft(int64_t N, int64_t k, int nb_bits);



typedef struct glwe_public_key {
  uint64_t N;
  uint64_t k;
  uint64_t Y;
  GLWECiphertext** pk;  // vector of Y element (A1,...,Ak, B)
  void* data;
} GLWEPublicKey;



// generation secret key
void glwe_secret_key_gen(const Core* core, uint64_t lambda, GLWESecretKey* s);

// generation public key
void glwe_public_key_gen(const Core* core, uint64_t lambda, GLWESecretKey* s,
                         GLWEPublicKey* pk);

#endif  // GLWE_KEY_H