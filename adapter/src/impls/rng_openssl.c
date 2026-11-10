#include "impls/rng_openssl.h"

#include <limits.h>
#include <openssl/evp.h>
#include <openssl/rand.h>
#include <stdint.h>
#include <string.h>

#include "logger.h"

// =============================================
// |                                           |
// |              System mode                  |
// |                                           |
// =============================================

static int openssl_system_bytes(const PvdaBackend* module, void* buf, size_t n)
{
	(void)module;
	uint8_t* p = buf;
	while (n > 0)
	{
		int len = n > INT_MAX ? INT_MAX : (int)n;
		if (RAND_bytes(p, len) != 1) return log_message(LOG_ERROR, "OpenSSL RAND_bytes failed");
		p += len;
		n -= (size_t)len;
	}
	return 0;
}

// =============================================
// |                                           |
// |          Seeded mode: AES-256-CTR         |
// |                                           |
// =============================================

// The keystream is the encryption of zeros: out = AES-CTR(0) = AES_key(IV), AES_key(IV + 1), ...
// The cipher context keeps the counter and the unused part of the current block between calls.
static int openssl_ctr_bytes(const PvdaBackend* module, void* buf, size_t n)
{
	EVP_CIPHER_CTX* ctx = module->rng_state;
	uint8_t* p          = buf;
	memset(p, 0, n);
	while (n > 0)
	{
		int len = n > (INT_MAX / 2) ? (INT_MAX / 2) : (int)n;
		int outl = 0;
		if (EVP_EncryptUpdate(ctx, p, &outl, p, len) != 1 || outl != len)
			return log_message(LOG_ERROR, "OpenSSL AES-256-CTR keystream generation failed");
		p += len;
		n -= (size_t)len;
	}
	return 0;
}

static void openssl_ctr_delete(void* state) { EVP_CIPHER_CTX_free(state); }

// =============================================
// |                                           |
// |              Init                         |
// |                                           |
// =============================================

int pvda_openssl_rng_init(PvdaRngVirtualTable* vt, void** state, const PvdaSeed* seed)
{
	if (!seed)
	{
		vt->pvda_rng_bytes = openssl_system_bytes;
		return 0;
	}

	EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
	if (!ctx) return log_message(LOG_ERROR, "OpenSSL cipher context allocation failed");

	// Seed conversion: the 32 seed bytes are the AES-256 key, the counter block starts at 0.
	static const uint8_t iv[16] = {0};
	if (EVP_EncryptInit_ex(ctx, EVP_aes_256_ctr(), NULL, seed->bytes, iv) != 1)
	{
		EVP_CIPHER_CTX_free(ctx);
		return log_message(LOG_ERROR, "OpenSSL AES-256-CTR initialization failed");
	}

	*state             = ctx;
	vt->pvda_rng_bytes = openssl_ctr_bytes;
	vt->delete_state   = openssl_ctr_delete;
	return 0;
}
