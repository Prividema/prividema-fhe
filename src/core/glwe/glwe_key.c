#include "glwe_key.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#include "glwe_params.h"
#include "logger.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "utils.h"

GLWESecretKey* alloc_glwe_secret_key(GLWEParams* params_glwe)
{
	uint64_t j        = 0;
	GLWESecretKey* sk = malloc(sizeof(GLWESecretKey));
	CHECK_ALLOC(sk, "sk's malloc failed in alloc_glwe_secret_key");

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

int uniform_glwe_secret_key(const MODULE* module, GLWESecretKey* sk, uint64_t nb_bits)
{
	uint64_t nn = pvda_module_extract_nn(module);
	// The Secret key values
	// Uniform random generation of k Zn[X] polynomials.
	CHECK_CALL(uniform_random_pol_znx(sk->values, nn * sk->k, nb_bits),
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

void delete_glwe_secret_key(GLWESecretKey* sk)
{
	if (!sk) return;
	free(sk->values);
	free(sk);
}

GLWESecretKeyDFT* alloc_glwe_secret_key_dft(GLWEParams* params_glwe)
{
	uint64_t j           = 0;
	GLWESecretKeyDFT* sk = malloc(sizeof(GLWESecretKeyDFT));
	CHECK_ALLOC(sk, "sk's malloc failed in new_glwe_secret_key.");
	sk->nn = params_glwe->nn;
	sk->k  = params_glwe->k;

	sk->values = malloc(params_glwe->nn * params_glwe->k * sizeof(PolyUnivDFT));
	CHECK_ALLOC(sk->values, "values' malloc failed in alloc_glwe_secret_key_dft");

	return sk;
cleanup:

	if (sk) free(sk->values);
	free(sk);
	return NULL;
}

PolyUnivDFT* glwe_sk_extract_poly_dft(const GLWESecretKeyDFT* sk_dft, uint64_t pos)
{
	assert(pos >= 0 && pos < sk_dft->k);
	return sk_dft->values + sk_dft->nn * pos;
}

void delete_glwe_secret_key_dft(GLWESecretKeyDFT* sk_dft)
{
	if (!sk_dft) return;

	free(sk_dft->values);
	free(sk_dft);
}
