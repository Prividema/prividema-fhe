#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <stdint.h>
#include <sys/types.h>

#include "backend.h"
#include "backend_rng.h"
#include "heint.h"

Test(heint_helpers, inverse_via_euclid)
{
	PvdaBackend* rng_backend = pvda_new_spqlios_backend(4);
	uint64_t primes[]        = {3, 5, 7, 13, 6701, 7793, 2147483647};

	for (size_t i = 0; i < sizeof(primes) / sizeof(primes[0]); ++i)
	{
		uint64_t pi = primes[i];
		for (size_t r = 0; r < 100; ++r)
		{
			uint64_t x;
			pvda_rand_uniform(rng_backend, &x, 1, pi - 1);

			uint64_t x_inv = mod_inv(x, pi);
			uint64_t prd   = (x * x_inv) % pi;

			cr_assert(eq(u64, prd, 1));
		}
	}

	pvda_delete_backend(rng_backend);
}

Test(heint_helpers, mont_mult)
{
	PvdaBackend* rng_backend = pvda_new_spqlios_backend(4);
	uint64_t primes[]        = {3, 5, 7, 13, 6701, 7793, 2147483647};

	for (size_t i = 0; i < sizeof(primes) / sizeof(primes[0]); ++i)
	{
		uint64_t pi     = primes[i];
		uint64_t p_tild = (1ull << 32) - mod_inv(pi, 1ull << 32);

		for (size_t r = 0; r < 100; ++r)
		{
			uint64_t x, y;
			pvda_rand_uniform(rng_backend, &x, 1, pi - 1);
			pvda_rand_uniform(rng_backend, &y, 1, pi - 1);

			uint64_t x_m, y_m;
			x_m = montgomery_encode_32bit(x, pi);
			y_m = montgomery_encode_32bit(y, pi);

			uint64_t xy_m        = montgomery_mult_32bit(x_m, y_m, pi, p_tild);
			uint64_t xy_observed = montgomery_decode_32bit(xy_m, pi, p_tild);

			uint64_t expc_val = (x * y) % pi;

			cr_assert(eq(u64, xy_observed, expc_val));
		}
	}

	pvda_delete_backend(rng_backend);
}
