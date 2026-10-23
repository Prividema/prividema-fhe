#include "glwe_key.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "backend.h"
#include "backend_arithmetic.h"
#include "glwe_params.h"
#include "rng.h"
#include "utils.h"

// =============================================
// |                                           |
// |                Secret Key                 |
// |                                           |
// =============================================

GLWESecretKey* new_glwe_secret_key(GLWEParams* params_glwe)
{
	uint64_t j        = 0;
	GLWESecretKey* sk = malloc(sizeof(GLWESecretKey));
	CHECK_ALLOC(sk, "sk's malloc failed in new_glwe_secret_key");

	sk->nn = params_glwe->nn;
	sk->k  = params_glwe->k;

	sk->values = calloc(params_glwe->nn * params_glwe->k, sizeof(PolyUniv));
	CHECK_ALLOC(sk->values, "values creation failed in glwe key generation");

	return sk;

cleanup:
	if (sk) free(sk->values);
	free(sk);
	return NULL;
}

void delete_glwe_secret_key(GLWESecretKey* sk)
{
	if (!sk) return;
	free(sk->values);
	free(sk);
}


int glwe_sk_random_uniform(const PvdaBackend* module, GLWESecretKey* sk, uint64_t nb_bits)
{
	uint64_t nn = pvda_module_extract_nn(module);
	// The Secret key values
	// Uniform random generation of k Zn[X] polynomials.
	CHECK_CALL(uniform_pow2_random_pol_znx(module, sk->values, nn * sk->k, nb_bits),
	           "random vector generation failed in key generation");

	return 0;

cleanup:
	return -1;
}

int glwe_sk_random_binary(const PvdaBackend* module, GLWESecretKey* sk)
{
	uint64_t nn = pvda_module_extract_nn(module);
	// The Secret key values
	// Uniform random generation of k Zn[X] polynomials.
	CHECK_CALL(binary_random_pol_znx(module, sk->values, nn * sk->k),
	           "random vector generation failed in key generation");

	return 0;

cleanup:
	return -1;
}

int glwe_sk_random_ternary(const PvdaBackend* module, GLWESecretKey* sk)
{
	uint64_t nn = pvda_module_extract_nn(module);
	// The Secret key values
	// Uniform random generation of k Zn[X] polynomials.
	CHECK_CALL(uniform_random_pol_znx(module, sk->values, nn * sk->k, -1, 1),
	           "random vector generation failed in key generation");

	return 0;

cleanup:
	return -1;
}

PolyUniv* glwe_sk_extract_poly(GLWESecretKey* sk, uint64_t pos)
{
	assert(pos >= 0 && pos < sk->k);
	return sk->values + sk->nn * pos;
}

// =============================================
// |                                           |
// |           Secret Key Prepared             |
// |                                           |
// =============================================

GLWESecretKeyPrepared* new_glwe_secret_key_prepared(const PvdaBackend* module, GLWEParams* params_glwe)
{
	GLWESecretKeyPrepared* sk = calloc(1, sizeof(GLWESecretKeyPrepared));
	CHECK_ALLOC(sk, "sk's allocation failed in new_glwe_secret_key_prepared");

	// Set before any other allocation: delete_glwe_secret_key_prepared needs it to free the components
	sk->module = module;
	sk->nn     = params_glwe->nn;
	sk->k      = params_glwe->k;

	// Coefficient-space copy of the key
	sk->values_coef = malloc(sk->nn * sk->k * sizeof(PolyUniv));
	CHECK_ALLOC(sk->values_coef, "values_coef's allocation failed in new_glwe_secret_key_prepared");

	// Array of k handles, zeroed so that a partial failure leaves NULLs that the delete function can skip
	sk->values_prep = calloc(sk->k, sizeof(BE_PolySvpPrep*));
	CHECK_ALLOC(sk->values_prep, "values_prep's allocation failed in new_glwe_secret_key_prepared");

	// One backend object per key component
	for (uint64_t j = 0; j < sk->k; ++j)
	{
		sk->values_prep[j] = pvda_new_svp_ppol(module);
		CHECK_ALLOC(sk->values_prep[j], "prepared key component allocation failed in new_glwe_secret_key_prepared");
	}

	return sk;

cleanup:
	delete_glwe_secret_key_prepared(sk);
	return NULL;
}

void delete_glwe_secret_key_prepared(GLWESecretKeyPrepared* sk_prep)
{
	if (!sk_prep) return;

	if (sk_prep->values_prep)
		for (uint64_t j = 0; j < sk_prep->k; ++j) pvda_delete_svp_ppol(sk_prep->module, sk_prep->values_prep[j]);

	free((void*)sk_prep->values_prep);
	free(sk_prep->values_coef);
	free(sk_prep);
}

BE_PolySvpPrep* glwe_prepared_sk_extract_poly_prep(const GLWESecretKeyPrepared* sk_prep, uint64_t pos)
{
	assert(pos < sk_prep->k);
	return sk_prep->values_prep[pos];
}

PolyUniv* glwe_sk_prepared_extract_poly_coefs(const GLWESecretKeyPrepared* sk_prep, uint64_t pos)
{
	assert(pos >= 0 && pos < sk_prep->k);
	return sk_prep->values_coef + sk_prep->nn * pos;
}

// =============================================
// |                                           |
// |           Secret Key Transforms           |
// |                                           |
// =============================================

int glwe_sk_prepare(const PvdaBackend* module, GLWESecretKeyPrepared* sk_prep, const GLWESecretKey* sk)
{
	int status  = -1;
	uint64_t nn = sk->nn;

	// The prepared key must have been allocated for the same key shape and with the same backend
	assert(sk_prep->nn == sk->nn && sk_prep->k == sk->k);
	assert(sk_prep->module == module);

	// Coefficient-space copy
	memcpy(sk_prep->values_coef, sk->values, nn * sk->k * sizeof(PolyUniv));

	// Each component is prepared as the scalar operand of scalar-vector products
	for (uint64_t j = 0; j < sk->k; ++j)
		CHECK_CALL(pvda_svp_prepare(module, sk_prep->values_prep[j], sk->values + j * nn),
		           "SVP preparation of a secret key component failed in glwe_sk_prepare");

	status = 0;

cleanup:
	return status;
}
