

#include "rng_openssl.h"

#include <assert.h>
#include <openssl/rand.h>
#include <stdint.h>

#include "backend.h"
#include "utils.h"

static inline void reduce_uniform_n(int64_t* tgt, int n_bits)
{
	int shft = 64 - n_bits;
	*tgt     = (int64_t)((uint64_t)(*tgt) << shft) >> shft;
}

static inline int openssl_read_rand(uint64_t* result, size_t bytes)
{
	size_t rand_bytes = RAND_bytes((uint8_t*)result, (int)bytes);
	if (rand_bytes != 1) return -1;
	return 0;
}
int openssl_rand_uniform_pow2(const PvdaBackend* module, int64_t* result, uint64_t nb_bits)
{
	// As result points to an uint64_t  nb_bits shall not exceed its size
	assert(nb_bits <= 8 * sizeof(int64_t));

	// If nb_bits equals the max. size, we just have to convert r to an int64_t.
	if (nb_bits == 8 * sizeof(int64_t))
		return openssl_read_rand((uint64_t*)result, 8);

	else
	{
		if (openssl_read_rand((uint64_t*)result, INT_ROUND_UP_DIV(nb_bits, 8)) < 0) return -1;

		reduce_uniform_n(result, (int)nb_bits);

		return 1;
	}

	return 0;
}

int openssl_rand_uniform(const PvdaBackend* module, int64_t* result, int64_t limit_down, int64_t limit_up)
{
	uint64_t max_delta = (uint64_t)limit_up - (uint64_t)limit_down;
	if (max_delta == UINT64_MAX) return -1;
	uint64_t bits = next_pow2_log(max_delta + 1);
	uint64_t mask = bits == 64 ? (UINT64_MAX) : (1ull << bits) - 1;

	uint64_t tmp = 0;

	int st;
	do
	{
		st = openssl_read_rand(&tmp, INT_ROUND_UP_DIV(bits, 8));
		if (st < 0) return -1;
		tmp &= mask;
	} while (tmp > max_delta);
	*result = (int64_t)((uint64_t)limit_down + tmp);
	return 0;
}

int openssl_rand_uniform_pow2_vec(const PvdaBackend* module, int64_t* res, uint64_t n, uint64_t nb_bits)
{
	CHECK_CALL(openssl_read_rand((uint64_t*)res, sizeof(int64_t) * n), "rng error");
	for (uint64_t p = 0; p < n; p++)
	{
		reduce_uniform_n(res + p, (int)nb_bits);
	}
	return 0;
cleanup:
	return -1;
}

int openssl_rand_uniform_binary_vec(const PvdaBackend* module, uint64_t* res, uint64_t n)
{
	CHECK_CALL(openssl_read_rand(res, sizeof(uint64_t) * n), "rng error");
	for (uint64_t p = 0; p < n; p++)
	{
		res[p] &= 1;
	}
	return 0;
cleanup:
	return -1;
}

void pvda_fill_openssl_rng(struct pvda_virtual_table* vt)
{
	vt->pvda_rand_uniform_pow2       = openssl_rand_uniform_pow2;
	vt->pvda_rand_uniform_pow2_vec   = openssl_rand_uniform_pow2_vec;
	vt->pvda_rand_uniform_binary_vec = openssl_rand_uniform_binary_vec;
	vt->pvda_rand_uniform            = openssl_rand_uniform;
}
