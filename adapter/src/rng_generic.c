#include <math.h>
#include <stdint.h>
#include <string.h>

#include "logger.h"
#include "backend_private.h"

/**
 * @file rng_generic.c
 *
 * Distributions derived from the random bytes of any RNG implementation. They are used by every implementation
 * that does not provide its own version of a distribution.
 */

// On some distros math.h doesn't define M_PI so we define it here just in case.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Number of 64-bit words drawn at once from the backend by the samplers below.
#define RNG_CHUNK_WORDS 512

// =============================================
// |                                           |
// |              Bit stream                   |
// |                                           |
// =============================================

/** @brief Buffered stream of random 64-bit words, refilled RNG_CHUNK_WORDS at a time. */
typedef struct
{
	const PvdaBackend* module;
	uint64_t buf[RNG_CHUNK_WORDS];
	size_t pos;
	size_t len;
} WordStream;

static inline void word_stream_init(WordStream* s, const PvdaBackend* module)
{
	s->module = module;
	s->pos = 0;
	s->len = 0;
}

/**
 * @brief Next word of the stream. remaining_hint is the number of words the caller still expects to need: refills
 * draw at most that many words (and at most RNG_CHUNK_WORDS), so small requests do not consume a whole chunk.
 */
static inline int word_stream_next(WordStream* s, uint64_t remaining_hint, uint64_t* out)
{
	if (s->pos == s->len)
	{
		size_t len = remaining_hint < RNG_CHUNK_WORDS ? (size_t)remaining_hint : RNG_CHUNK_WORDS;
		if (len == 0) len = 1;
		if (s->module->rng_vt.pvda_rng_bytes(s->module, s->buf, len * sizeof(uint64_t)) < 0) return -1;
		s->pos = 0;
		s->len = len;
	}
	*out = s->buf[s->pos++];
	return 0;
}

// =============================================
// |                                           |
// |              Uniform                      |
// |                                           |
// =============================================

static int generic_uniform_pow2_vec(const PvdaBackend* module, int64_t* res, uint64_t n, uint64_t nb_bits)
{
	if (nb_bits == 0 || nb_bits > 64)
	{
		log_message(LOG_ERROR, "uniform_pow2: nb_bits must be in [1, 64]");
		return -1;
	}
	if (n == 0) return 0;
	if (module->rng_vt.pvda_rng_bytes(module, res, n * sizeof(int64_t)) < 0) return -1;

	if (nb_bits < 64)
	{
		// Keep the nb_bits low bits, sign-extended: uniform in [-2^(nb_bits-1), 2^(nb_bits-1))
		int shift = 64 - (int)nb_bits;
		for (uint64_t i = 0; i < n; ++i) res[i] = (int64_t)((uint64_t)res[i] << shift) >> shift;
	}
	return 0;
}

static int generic_uniform_range_vec(const PvdaBackend* module, int64_t* res, uint64_t n, int64_t low, int64_t high)
{
	if (low > high)
	{
		log_message(LOG_ERROR, "uniform_range: low must be <= high");
		return -1;
	}
	uint64_t delta = (uint64_t)high - (uint64_t)low;  // number of values - 1

	// Full range: every 64-bit word is a valid sample
	if (delta == UINT64_MAX) return module->rng_vt.pvda_rng_bytes(module, res, n * sizeof(int64_t));

	// Rejection sampling on the smallest power of 2 covering [0, delta]: acceptance probability > 1/2
	int bits      = delta == 0 ? 0 : 64 - __builtin_clzll(delta);
	uint64_t mask = bits == 64 ? UINT64_MAX : (1ULL << bits) - 1;

	WordStream s;
	word_stream_init(&s, module);
	for (uint64_t i = 0; i < n; ++i)
	{
		uint64_t w;
		do
		{
			if (word_stream_next(&s, n - i, &w) < 0) return -1;
			w &= mask;
		} while (w > delta);
		res[i] = (int64_t)((uint64_t)low + w);
	}
	return 0;
}

static int generic_binary_vec(const PvdaBackend* module, int64_t* res, uint64_t n)
{
	WordStream s;
	word_stream_init(&s, module);
	uint64_t w = 0;
	for (uint64_t i = 0; i < n; ++i)
	{
		if (i % 64 == 0 && word_stream_next(&s, (n - i + 63) / 64, &w) < 0) return -1;
		res[i] = (int64_t)(w & 1);
		w >>= 1;
	}
	return 0;
}

// =============================================
// |                                           |
// |              Binomial                     |
// |                                           |
// =============================================

/** @brief Number of ones among k fresh random bits of the stream. */
static inline int popcount_bits(WordStream* s, uint64_t k, uint64_t words_hint, uint64_t* out)
{
	uint64_t acc = 0;
	uint64_t w;
	while (k >= 64)
	{
		if (word_stream_next(s, words_hint, &w) < 0) return -1;
		acc += (uint64_t)__builtin_popcountll(w);
		k -= 64;
	}
	if (k)
	{
		if (word_stream_next(s, words_hint, &w) < 0) return -1;
		acc += (uint64_t)__builtin_popcountll(w & ((1ULL << k) - 1));
	}
	*out = acc;
	return 0;
}

static int generic_centered_binomial_vec(const PvdaBackend* module, int64_t* res, uint64_t n, uint64_t eta)
{
	if (eta > (uint64_t)INT64_MAX)
	{
		log_message(LOG_ERROR, "centered_binomial: eta too large");
		return -1;
	}
	if (eta == 0)
	{
		memset(res, 0, n * sizeof(int64_t));
		return 0;
	}

	uint64_t words_per_sample = 2 * ((eta + 63) / 64);
	WordStream s;
	word_stream_init(&s, module);
	for (uint64_t i = 0; i < n; ++i)
	{
		uint64_t hint = (n - i) * words_per_sample;
		uint64_t a, b;
		if (popcount_bits(&s, eta, hint, &a) < 0) return -1;
		if (popcount_bits(&s, eta, hint, &b) < 0) return -1;
		res[i] = (int64_t)a - (int64_t)b;
	}
	return 0;
}

// =============================================
// |                                           |
// |              Normal                       |
// |                                           |
// =============================================

static int generic_normal_vec(const PvdaBackend* module, double* res, uint64_t n, double mu, double sigma)
{
	if (!isfinite(mu) || !isfinite(sigma) || sigma < 0)
	{
		log_message(LOG_ERROR, "normal: mu and sigma must be finite and sigma >= 0");
		return -1;
	}

	// Box-Muller: two uniforms give two independent standard normal samples.
	WordStream s;
	word_stream_init(&s, module);
	for (uint64_t i = 0; i < n; i += 2)
	{
		uint64_t w1, w2;
		uint64_t hint = n - i + 1;
		if (word_stream_next(&s, hint, &w1) < 0) return -1;
		if (word_stream_next(&s, hint, &w2) < 0) return -1;

		double u1 = (double)((w1 >> 11) + 1) * 0x1p-53;  // (0, 1]: log(u1) is finite
		double u2 = (double)(w2 >> 11) * 0x1p-53;        // [0, 1)

		double r     = sqrt(-2.0 * log(u1));
		double theta = 2.0 * M_PI * u2;

		res[i] = mu + sigma * r * cos(theta);
		if (i + 1 < n) res[i + 1] = mu + sigma * r * sin(theta);
	}
	return 0;
}

void pvda_rng_fill_generic(PvdaRngVirtualTable* vt)
{
	vt->pvda_rng_uniform_pow2_vec      = generic_uniform_pow2_vec;
	vt->pvda_rng_uniform_range_vec     = generic_uniform_range_vec;
	vt->pvda_rng_binary_vec            = generic_binary_vec;
	vt->pvda_rng_centered_binomial_vec = generic_centered_binomial_vec;
	vt->pvda_rng_normal_vec            = generic_normal_vec;
}
