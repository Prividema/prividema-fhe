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

typedef struct glwe_secret_key {
  uint64_t n;  // n= N*k (n is constant: N and k can change when we change the
               // module N)
  int64_t* values;
  void* data;
} GLWESecretKey;

GLWESecretKey* new_glwe_secret_key(int64_t* values, int64_t n);
void delete_glwe_secret_key(GLWESecretKey* sk);

typedef struct glwe_prep_secret_key {
  uint64_t N;
  uint64_t k;
  PolyUnivDFT** values;  // vec of size k, each element is prepared vec
  void* data;
} GLWEPreparedSK;

/**
 * @brief Creates a GLWE Secret key in DFT space
 * 
 * @param values The values of the secret key in DFT space.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of ZnX polynomial in the secret key.
 * @return GLWEPreparedSK* 
 */
GLWEPreparedSK* new_glwe_secret_key_dft(PolyUnivDFT** values, int64_t N, int64_t k);

/**
 * @brief Delete the secret key that is in DFT space.
 * 
 * @param sk_dft The secret key in DFT space.
 */
void delete_glwe_secret_key_dft(GLWEPreparedSK* sk_dft);

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