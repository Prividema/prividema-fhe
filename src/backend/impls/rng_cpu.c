
#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <bcrypt.h>
#include <windows.h>
#pragma comment(lib, "bcrypt.lib")
#endif
#ifdef __linux__
#include <sys/random.h>
#endif

#include <assert.h>
#include <stdint.h>

#include "backend.h"
#include "backend_private.h"
#include "utils.h"

/*
    Read a random number depending on the OS :
    - On Windows : Uses Windows' Cryptographic API called CNG.
    - On MACOS/FreeBSD : Call to arc4random_buf.
                        According to arc4random's doc, the whole program crashes
                        if an error occurs during the generation.
    - On other Linux distributions : read /dev/urandom.
*/
static inline int cpu_read_rand(uint64_t* result, size_t bytes)
{
// For Windows
#ifdef _WIN32
	NTSTATUS status = BCryptGenRandom(NULL, (PUCHAR)result, bytes, BCRYPT_USE_SYSTEM_PREFERRED_RNG);
	if (status != STATUS_SUCCESS) return log_message(LOG_ERROR, "BCryptGenRandom() Failed");

// For MACOS/FreeBSD
// According to arc4random's doc, the function crashes if an error occurs :
// "Cryptographic randomness is considered fundamental — if it’s broken, continuing execution is unsafe."
#elif defined(__APPLE__) || defined(__FreeBSD__)
	arc4random_buf(result, bytes);

// For Linux
#elif defined(__linux__)

	size_t rand_bytes = getrandom(result, bytes, 0);
	if (rand_bytes != bytes) return -1;
	return 0;

// I don't know what system this block below would be (and it would be quite slow)
// But for now we leave it for compatibility
#else
	// THIS IS VERY SLOW!
	FILE* f = fopen("/dev/urandom", "rb");
	if (!f) return log_perror("fopen");
	int r = fread(result, 1, bytes, f);
	fclose(f);
	if (r != bytes) return log_perror("fread");
#endif

	return 0;
}

static inline void reduce_uniform_n(int64_t* tgt, int n_bits)
{
	int shft = 64 - n_bits;
	*tgt     = (int64_t)((uint64_t)(*tgt) << shft) >> shft;
}

/**
 * @brief Generates a uniformly sampled random number in [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param module  The backend object
 * @param result  The resulting uniformly sampled integer
 * @param nb_bits The number of bits of the result
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ref_rand_uniform_pow2(const PvdaBackend* module, int64_t* result, uint64_t nb_bits)
{
	// As result points to an uint64_t  nb_bits shall not exceed its size
	assert(nb_bits <= 8 * sizeof(int64_t));

	// If nb_bits equals the max. size, we just have to convert r to an int64_t.
	if (nb_bits == 8 * sizeof(int64_t))
		return cpu_read_rand((uint64_t*)result, 8);

	else
	{
		if (cpu_read_rand((uint64_t*)result, INT_ROUND_UP_DIV(nb_bits, 8)) < 0) return -1;

		reduce_uniform_n(result, (int)nb_bits);

		return 1;
	}

	return 0;
}

/**
 * @brief Generates a uniformly sampled random number in [limit_down, limit_up] via
 * power-of-2 sampling and resampling if out-of-bounds
 *
 * @param module      The backend object
 * @param result      The resulting uniformly sampled integer
 * @param limit_down  The lower bound of the uniform sample
 * @param limit_up    The upper bound of the uniform sample
 *
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ref_rand_uniform(const PvdaBackend* module, int64_t* result, int64_t limit_down, int64_t limit_up)
{
	uint64_t max_delta = (uint64_t)limit_up - (uint64_t)limit_down;
	if (max_delta == UINT64_MAX) return -1;
	uint64_t bits = next_pow2_log(max_delta + 1);
	uint64_t mask = bits == 64 ? (UINT64_MAX) : (1ull << bits) - 1;

	uint64_t tmp = 0;

	int st;
	do
	{
		st = cpu_read_rand(&tmp, INT_ROUND_UP_DIV(bits, 8));
		if (st < 0) return -1;
		tmp &= mask;
	} while (tmp > max_delta);
	*result = (int64_t)((uint64_t)limit_down + tmp);
	return 0;
}

/**
 * @brief Generates a uniformly random vector of numbers
 *
 * Coefficients are uniformly sampled in range [-2^(nb_bits-1), 2^(nb_bits-1))
 *
 * @param module  The backend object
 * @param res     The result uniformly drawn vector of numbers
 * @param nn      Number of coeffients in the polynomial (eq. degree of the cyclotomial poly)
 * @param nb_bits Number of randomness bits per coefficient.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int ref_rand_uniform_pow2_vec(const PvdaBackend* module, int64_t* res, uint64_t n, uint64_t nb_bits)
{
	CHECK_CALL(cpu_read_rand((uint64_t*)res, sizeof(int64_t) * n), "rng error");
	for (uint64_t p = 0; p < n; p++)
	{
		reduce_uniform_n(res + p, (int)nb_bits);
	}
	return 0;
cleanup:
	return -1;
}

/**
 * @brief Generates a uniformly random binary vector of numbers
 *
 * Coefficients are uniformly sampled in {0, 1}
 *
 * @param module  The backend object
 * @param res     The result uniformly drawn \ZnX polynomial.
 * @param nn      Number of coeffients in the polynomial (eq. degree of the cyclotomial poly)
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int ref_rand_uniform_binary_vec(const PvdaBackend* module, uint64_t* res, uint64_t n)
{
	CHECK_CALL(cpu_read_rand(res, sizeof(uint64_t) * n), "rng error");
	for (uint64_t p = 0; p < n; p++)
	{
		res[p] &= 1;
	}
	return 0;
cleanup:
	return -1;
}

void pvda_fill_ref_rng(struct pvda_virtual_table* vt)
{
	vt->pvda_rand_uniform_pow2       = ref_rand_uniform_pow2;
	vt->pvda_rand_uniform_pow2_vec   = ref_rand_uniform_pow2_vec;
	vt->pvda_rand_uniform_binary_vec = ref_rand_uniform_binary_vec;
	vt->pvda_rand_uniform            = ref_rand_uniform;
}
