
#include <criterion/criterion.h>
#include <criterion/internal/new_asserts.h>
#include <criterion/new/assert.h>
#include <stdint.h>

#include "backend.h"
#include "backend_arithmetic.h"
#include "backend_nt.h"
#include "backend_rng.h"
#include "core/glwe/glwe_params.h"
#include "rng.h"
#include "stat_utils.h"
#include "test_utils.h"

static inline uint64_t pow_exp_mod(uint64_t base, uint64_t exp, uint64_t m)
{
	uint64_t ans  = 1;
	uint64_t mult = base;

	while (exp > 0)
	{
		if (exp & 1)
		{
			ans *= mult;
			ans %= m;
		}
		exp  = exp / 2;
		mult = (mult * mult) % m;
	}

	return ans;
}

/**
 * Tetsing that the algorithm to find the 2N'th negative root of unity actually finds a
 * 2N'th negative root of unity
 */
Test(twoth_root, works_16)
{
	uint64_t nn          = 16;
	PvdaBackend* backend = pvda_new_spqlios_backend(nn);

	uint64_t ps[] = {3329, 12289, 8380417};

	for (int i = 0; i < sizeof(ps) / sizeof(ps[0]); ++i)
	{
		uint64_t root = pvda_ring_2nth_root(backend, ps[i]);
		uint64_t res  = pow_exp_mod(root, nn, ps[i]);
		cr_assert(eq(i64, res, ps[i] - 1ll), "Error with root %d", root);
	}
	// cr_assert()
	//
	pvda_delete_backend(backend);
}
