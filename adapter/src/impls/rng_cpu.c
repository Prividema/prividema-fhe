#include "impls/rng_cpu.h"

#ifdef _WIN32
#include <bcrypt.h>
#include <windows.h>
#pragma comment(lib, "bcrypt.lib")
#endif
#ifdef __linux__
#include <errno.h>
#include <sys/random.h>
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "logger.h"

// =============================================
// |                                           |
// |              System mode                  |
// |                                           |
// =============================================

/*
    Read random bytes depending on the OS :
    - On Windows : Uses Windows' Cryptographic API called CNG.
    - On MACOS/FreeBSD : Call to arc4random_buf.
                        According to arc4random's doc, the whole program crashes
                        if an error occurs during the generation.
    - On Linux : Call to getrandom().
    - On other systems : read /dev/urandom (slow fallback).
*/
int cpu_os_random_bytes(void* buf, size_t n)
{
#ifdef _WIN32
	NTSTATUS status = BCryptGenRandom(NULL, (PUCHAR)buf, (ULONG)n, BCRYPT_USE_SYSTEM_PREFERRED_RNG);
	if (status != STATUS_SUCCESS) return log_message(LOG_ERROR, "BCryptGenRandom() Failed");
	return 0;

#elif defined(__APPLE__) || defined(__FreeBSD__)
	// "Cryptographic randomness is considered fundamental — if it’s broken, continuing execution is unsafe."
	arc4random_buf(buf, n);
	return 0;

#elif defined(__linux__)
	// getrandom can return fewer bytes than requested (large requests, signals): loop until done.
	uint8_t* p = buf;
	while (n > 0)
	{
		ssize_t r = getrandom(p, n, 0);
		if (r < 0)
		{
			if (errno == EINTR) continue;
			log_perror("getrandom");
			return -1;
		}
		p += r;
		n -= (size_t)r;
	}
	return 0;

#else
	// Fallback for unknown systems
	// THIS IS VERY SLOW!
	FILE* f = fopen("/dev/urandom", "rb");
	if (!f) return log_perror("fopen");
	size_t r = fread(buf, 1, n, f);
	fclose(f);
	if (r != n) return log_perror("fread");
	return 0;
#endif
}

static int cpu_system_bytes(const PvdaBackend* module, void* buf, size_t n)
{
	(void)module;
	return cpu_os_random_bytes(buf, n);
}

// =============================================
// |                                           |
// |              Seeded mode: ChaCha20        |
// |                                           |
// =============================================

static inline uint32_t rotl32(uint32_t x, int r) { return (x << r) | (x >> (32 - r)); }

#define CHACHA_QR(a, b, c, d)  \
	do                         \
	{                          \
		a += b;                \
		d = rotl32(d ^ a, 16); \
		c += d;                \
		b = rotl32(b ^ c, 12); \
		a += b;                \
		d = rotl32(d ^ a, 8);  \
		c += d;                \
		b = rotl32(b ^ c, 7);  \
	} while (0)

void pvda_chacha20_block(const uint32_t key[8], uint64_t counter, uint64_t nonce, uint8_t out[64])
{
	// "expand 32-byte k"
	uint32_t in[16] = {0x61707865,
	                   0x3320646e,
	                   0x79622d32,
	                   0x6b206574,
	                   key[0],
	                   key[1],
	                   key[2],
	                   key[3],
	                   key[4],
	                   key[5],
	                   key[6],
	                   key[7],
	                   (uint32_t)counter,
	                   (uint32_t)(counter >> 32),
	                   (uint32_t)nonce,
	                   (uint32_t)(nonce >> 32)};
	uint32_t x[16];
	memcpy(x, in, sizeof(x));

	for (int i = 0; i < 10; ++i)
	{
		// Column rounds
		CHACHA_QR(x[0], x[4], x[8], x[12]);
		CHACHA_QR(x[1], x[5], x[9], x[13]);
		CHACHA_QR(x[2], x[6], x[10], x[14]);
		CHACHA_QR(x[3], x[7], x[11], x[15]);
		// Diagonal rounds
		CHACHA_QR(x[0], x[5], x[10], x[15]);
		CHACHA_QR(x[1], x[6], x[11], x[12]);
		CHACHA_QR(x[2], x[7], x[8], x[13]);
		CHACHA_QR(x[3], x[4], x[9], x[14]);
	}

	for (int i = 0; i < 16; ++i)
	{
		uint32_t v     = x[i] + in[i];
		out[4 * i + 0] = (uint8_t)v;
		out[4 * i + 1] = (uint8_t)(v >> 8);
		out[4 * i + 2] = (uint8_t)(v >> 16);
		out[4 * i + 3] = (uint8_t)(v >> 24);
	}
}

typedef struct
{
	uint32_t key[8];
	uint64_t counter;   ///< Next block to generate
	uint8_t block[64];  ///< Current keystream block
	size_t pos;         ///< Bytes of block already used
} ChaChaState;

static int cpu_chacha_bytes(const PvdaBackend* module, void* buf, size_t n)
{
	ChaChaState* st = module->rng_state;
	uint8_t* out    = buf;
	while (n > 0)
	{
		if (st->pos == 64)
		{
			pvda_chacha20_block(st->key, st->counter++, 0, st->block);
			st->pos = 0;
		}
		size_t take = 64 - st->pos < n ? 64 - st->pos : n;
		memcpy(out, st->block + st->pos, take);
		st->pos += take;
		out += take;
		n -= take;
	}
	return 0;
}

static void cpu_chacha_delete(void* state)
{
	if (!state) return;
	// Do not leave the key in freed memory
	volatile uint8_t* p = state;
	for (size_t i = 0; i < sizeof(ChaChaState); ++i) p[i] = 0;
	free(state);
}

// =============================================
// |                                           |
// |              Init                         |
// |                                           |
// =============================================

int pvda_cpu_rng_init(PvdaRngVirtualTable* vt, void** state, const PvdaSeed* seed)
{
	if (!seed)
	{
		vt->pvda_rng_bytes = cpu_system_bytes;
		return 0;
	}

	ChaChaState* st = calloc(1, sizeof(ChaChaState));
	if (!st) return log_perror("ChaCha20 state allocation failed");

	// Seed conversion: the 32 seed bytes are the ChaCha20 key, as 8 little-endian words.
	for (int i = 0; i < 8; ++i)
	{
		const uint8_t* b = seed->bytes + 4 * i;
		st->key[i]       = (uint32_t)b[0] | (uint32_t)b[1] << 8 | (uint32_t)b[2] << 16 | (uint32_t)b[3] << 24;
	}
	st->counter = 0;
	st->pos     = 64;  // no block generated yet

	*state             = st;
	vt->pvda_rng_bytes = cpu_chacha_bytes;
	vt->delete_state   = cpu_chacha_delete;
	return 0;
}
