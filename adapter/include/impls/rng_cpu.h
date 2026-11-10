#ifndef RNG_CPU_H
#define RNG_CPU_H

#include <stddef.h>

#include "backend_private.h"

/**
 * @file rng_cpu.h
 *
 * CPU RNG implementation: system randomness of the OS, or ChaCha20 when seeded.
 */

/**
 * @brief Reads n bytes from the randomness source of the OS
 * (getrandom on Linux, arc4random_buf on macOS/FreeBSD, BCryptGenRandom on Windows, /dev/urandom otherwise).
 */
int cpu_os_random_bytes(void* buf, size_t n);

/**
 * @brief Initializes an RNG virtual table and state as the CPU RNG implementation.
 *
 * @param vt    The virtual table, which already holds the generic distributions.
 * @param state Receives the RNG state (left NULL in system mode).
 * @param seed  NULL for system mode, the seed for seeded (ChaCha20) mode.
 *
 * The seed is converted into a ChaCha20 key: the 32 seed bytes, read as eight little-endian 32-bit words.
 * The nonce is 0 and the 64-bit block counter starts at 0 (original ChaCha20 layout).
 */
int pvda_cpu_rng_init(PvdaRngVirtualTable* vt, void** state, const PvdaSeed* seed);

/**
 * @brief Computes one 64-byte ChaCha20 block (20 rounds) for the given key, 64-bit counter and 64-bit nonce.
 *
 * Exposed for testing against reference vectors.
 */
void pvda_chacha20_block(const uint32_t key[8], uint64_t counter, uint64_t nonce, uint8_t out[64]);

#endif
