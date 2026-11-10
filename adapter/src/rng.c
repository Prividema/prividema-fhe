#include <stdlib.h>

#include "backend_private.h"
#include "backend_rng.h"
#include "impls/rng_cpu.h"
#include "logger.h"

#ifdef OPENSSL_BACKEND
#include "impls/rng_openssl.h"
#endif

// =============================================
// |                                           |
// |              Creation                     |
// |                                           |
// =============================================

int pvda_rng_type_available(PvdaRngType type)
{
	switch (type)
	{
		case PVDA_RNG_CPU:
			return 1;
		case PVDA_RNG_OPENSSL:
#ifdef OPENSSL_BACKEND
			return 1;
#else
			return 0;
#endif
	}
	return 0;
}

PvdaRngType pvda_default_rng_type(void)
{
	return pvda_rng_type_available(PVDA_RNG_OPENSSL) ? PVDA_RNG_OPENSSL : PVDA_RNG_CPU;
}

int pvda_backend_set_rng(PvdaBackend* module, PvdaRngType type, const PvdaSeed* seed)
{
	if (!pvda_rng_type_available(type)) return log_message(LOG_ERROR, "Requested RNG is not available in this build");

	// Build the new RNG aside, so that the current one is kept if the initialization fails
	PvdaRngVirtualTable vt = {0};
	void* state            = NULL;
	pvda_rng_fill_generic(&vt);

	int status = -1;
	switch (type)
	{
		case PVDA_RNG_CPU:
			status = pvda_cpu_rng_init(&vt, &state, seed);
			break;
		case PVDA_RNG_OPENSSL:
#ifdef OPENSSL_BACKEND
			status = pvda_openssl_rng_init(&vt, &state, seed);
#endif
			break;
	}

	if (status < 0 || !vt.pvda_rng_bytes)
	{
		if (vt.delete_state) vt.delete_state(state);
		return log_message(LOG_ERROR, "RNG initialization failed");
	}

	pvda_delete_backend_rng(module);
	module->rng_vt     = vt;
	module->rng_type   = type;
	module->rng_seeded = seed != NULL;
	module->rng_state  = state;
	return 0;
}

void pvda_delete_backend_rng(PvdaBackend* module)
{
	if (module->rng_vt.delete_state) module->rng_vt.delete_state(module->rng_state);
	module->rng_vt    = (PvdaRngVirtualTable){0};
	module->rng_state = NULL;
}

PvdaRngType pvda_backend_rng_type(const PvdaBackend* module) { return module->rng_type; }

int pvda_backend_rng_is_seeded(const PvdaBackend* module) { return module->rng_seeded; }

int pvda_seed_generate(PvdaSeed* seed) { return cpu_os_random_bytes(seed->bytes, PVDA_SEED_BYTES); }

// =============================================
// |                                           |
// |              Distributions                |
// |                                           |
// =============================================

int pvda_rng_bytes(const PvdaBackend* module, void* buf, size_t n) { return module->rng_vt.pvda_rng_bytes(module, buf, n); }

int pvda_rng_uniform_pow2_vec(const PvdaBackend* module, int64_t* res, uint64_t n, uint64_t nb_bits)
{
	return module->rng_vt.pvda_rng_uniform_pow2_vec(module, res, n, nb_bits);
}

int pvda_rng_uniform_range_vec(const PvdaBackend* module, int64_t* res, uint64_t n, int64_t low, int64_t high)
{
	return module->rng_vt.pvda_rng_uniform_range_vec(module, res, n, low, high);
}

int pvda_rng_binary_vec(const PvdaBackend* module, int64_t* res, uint64_t n) { return module->rng_vt.pvda_rng_binary_vec(module, res, n); }

int pvda_rng_centered_binomial_vec(const PvdaBackend* module, int64_t* res, uint64_t n, uint64_t eta)
{
	return module->rng_vt.pvda_rng_centered_binomial_vec(module, res, n, eta);
}

int pvda_rng_normal_vec(const PvdaBackend* module, double* res, uint64_t n, double mu, double sigma)
{
	return module->rng_vt.pvda_rng_normal_vec(module, res, n, mu, sigma);
}
