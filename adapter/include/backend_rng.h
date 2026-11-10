#ifndef PVDA_BACK_RNG_H
#define PVDA_BACK_RNG_H

#include <stddef.h>
#include <stdint.h>

#include "backend.h"
#include "pvda_types.h"

/**
 * @file backend_rng.h
 *
 * @brief Random number generation of a backend.
 *
 * Every PvdaBackend holds, next to its arithmetic virtual table, an RNG virtual table and its state: all the
 * randomness used by the library with a backend (key generation, encryption masks and noise, ...) comes from it.
 * The RNG implementation is chosen independently of the arithmetic implementation (pvda_backend_set_rng).
 *
 * Each RNG implementation can run in two modes:
 * - system mode (seed == NULL): non-deterministic, the randomness comes from the system source of the
 *   implementation;
 * - seeded mode: deterministic, the same seed always gives the same random stream.
 *
 * | RNG type          | System mode                                          | Seeded mode   |
 * |-------------------|------------------------------------------------------|---------------|
 * | PVDA_RNG_CPU      | OS source (getrandom, arc4random, BCryptGenRandom)   | ChaCha20      |
 * | PVDA_RNG_OPENSSL  | OpenSSL RAND_bytes                                   | AES-256-CTR   |
 *
 * PVDA_RNG_OPENSSL is only available when the project is configured with -DOPENSSL_BACKEND=ON.
 * A new backend uses pvda_default_rng_type() in system mode.
 *
 * Every RNG type provides all the distributions below. An implementation only has to provide random bytes: the
 * distributions are then derived from them by the adapter, unless the implementation has its own version.
 *
 * @warning The RNG state of a backend is not thread-safe: do not sample from the same backend in several threads.
 * @warning The samplers are not constant-time.
 *
 * All functions returning an int return 0 on success and -1 on failure (including invalid parameters).
 */

/** @brief The available RNG implementations */
typedef enum pvda_rng_type_t
{
	PVDA_RNG_CPU,      ///< System randomness of the OS, ChaCha20 when seeded
	PVDA_RNG_OPENSSL,  ///< OpenSSL RAND_bytes, AES-256-CTR when seeded
} PvdaRngType;

// =============================================
// |                                           |
// |              Configuration                |
// |                                           |
// =============================================

/**
 * @brief Tells whether an RNG implementation has been compiled in.
 *
 * @retval 1 if the implementation is available, 0 otherwise.
 */
int pvda_rng_type_available(PvdaRngType type);

/** @brief The RNG type of a new backend: PVDA_RNG_OPENSSL if available, PVDA_RNG_CPU otherwise. */
PvdaRngType pvda_default_rng_type(void);

/**
 * @brief Replaces the RNG of a backend.
 *
 * The seed is converted by the adapter into the seed format of the implementation. The same (type, seed) pair
 * always produces the same random stream, whatever the arithmetic implementation of the backend. Calling it again
 * with the same seed restarts the stream from the beginning.
 *
 * @param module The backend.
 * @param type   The RNG implementation.
 * @param seed   NULL for system mode (non-deterministic), the seed for seeded mode (deterministic).
 *
 * @retval -1 if type is not available or in case of failure; the previous RNG of the backend is then kept.
 */
int pvda_backend_set_rng(PvdaBackend* module, PvdaRngType type, const PvdaSeed* seed);

/** @brief Returns the RNG type of a backend. */
PvdaRngType pvda_backend_rng_type(const PvdaBackend* module);

/** @brief Tells whether the RNG of a backend is seeded (deterministic). @retval 1 if seeded, 0 otherwise. */
int pvda_backend_rng_is_seeded(const PvdaBackend* module);

/**
 * @brief Generates a fresh seed from the system randomness source (the library generates the seed).
 *
 * Save it to reproduce the random stream later with pvda_backend_set_rng.
 */
int pvda_seed_generate(PvdaSeed* seed);

// =============================================
// |                                           |
// |              Distributions                |
// |                                           |
// =============================================

/**
 * @brief Fills buf with n uniformly random bytes.
 *
 * The module is const as for the arithmetic functions, but its RNG state advances.
 */
int pvda_rng_bytes(const PvdaBackend* module, void* buf, size_t n);

/**
 * @brief Samples n integers uniformly in [-2^(nb_bits-1), 2^(nb_bits-1)).
 *
 * @pre 1 <= nb_bits <= 64 (returns -1 otherwise). With nb_bits = 64, all int64_t values are possible.
 */
int pvda_rng_uniform_pow2_vec(const PvdaBackend* module, int64_t* res, uint64_t n, uint64_t nb_bits);

/**
 * @brief Samples n integers uniformly in [low, high] (both included), without bias (rejection sampling).
 *
 * @pre low <= high (returns -1 otherwise). The full int64_t range is allowed.
 */
int pvda_rng_uniform_range_vec(const PvdaBackend* module, int64_t* res, uint64_t n, int64_t low, int64_t high);

/** @brief Samples n integers uniformly in {0, 1}. */
int pvda_rng_binary_vec(const PvdaBackend* module, int64_t* res, uint64_t n);

/**
 * @brief Samples n integers from the centered binomial distribution of parameter eta.
 *
 * Each sample is (b_1 + ... + b_eta) - (b'_1 + ... + b'_eta) for 2 * eta independent uniform bits, i.e.
 * B(2 * eta, 1/2) - eta. The samples lie in [-eta, eta], with mean 0 and variance eta / 2.
 * eta = 0 gives only zeros. The cost is O(eta) random bits per sample.
 */
int pvda_rng_centered_binomial_vec(const PvdaBackend* module, int64_t* res, uint64_t n, uint64_t eta);

/**
 * @brief Samples n reals from the normal distribution of mean mu and standard deviation sigma.
 *
 * Uses the Box-Muller transform on 53-bit uniform doubles.
 *
 * @pre sigma >= 0, mu and sigma finite (returns -1 otherwise).
 */
int pvda_rng_normal_vec(const PvdaBackend* module, double* res, uint64_t n, double mu, double sigma);

#endif
