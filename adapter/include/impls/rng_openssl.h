#ifndef RNG_OPENSSL_H
#define RNG_OPENSSL_H

#include "backend_private.h"

/**
 * @file rng_openssl.h
 *
 * OpenSSL RNG implementation: RAND_bytes, or AES-256-CTR when seeded.
 */

/**
 * @brief Initializes an RNG virtual table and state as the OpenSSL RNG implementation.
 *
 * @param vt    The virtual table, which already holds the generic distributions.
 * @param state Receives the RNG state (left NULL in system mode).
 * @param seed  NULL for system mode, the seed for seeded (AES-256-CTR) mode.
 *
 * The seed is converted into an AES-256 key (the 32 seed bytes). The 128-bit counter block (IV) starts at 0, and
 * the random stream is the AES-256-CTR keystream, i.e. AES_key(0), AES_key(1), ...
 */
int pvda_openssl_rng_init(PvdaRngVirtualTable* vt, void** state, const PvdaSeed* seed);

#endif
