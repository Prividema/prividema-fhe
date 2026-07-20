
#include <assert.h>
#include <stdint.h>

#include "backend.h"
#include "backend_private.h"
#include "backend_rng.h"

uint64_t pow_exp(uint64_t base, uint64_t exp, uint64_t m)
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

uint64_t ref_ring_2nth_root(const PvdaBackend* backend, uint64_t t)
{
	uint64_t nn = pvda_module_extract_nn(backend);

	assert(t <= INT32_MAX);
	assert((t - 1) % (2 * nn) == 0);  // 2n must divide t-1

	while (1)
	{
		uint64_t candidate;
		pvda_rand_uniform(backend, &candidate, 0, t - 1);  //TODO: check_call

		uint64_t g = pow_exp(candidate, (t - 1) / (2 * nn), t);
		if (pow_exp(g, nn, t) == t - 1)
		{
			return g;
		}
	}
}

void pvda_fill_ref_nt(struct pvda_virtual_table* vt) { vt->pvda_ring_2nth_root = ref_ring_2nth_root; }
