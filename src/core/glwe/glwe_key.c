#include "glwe_key.h"

#include <assert.h>
#include <stdint.h>

#include "logger.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "utils.h"

GLWESecretKey* alloc_glwe_secret_key(uint64_t N, uint64_t k)
{
	uint64_t j        = 0;
	GLWESecretKey* sk = malloc(sizeof(GLWESecretKey));
	CHECK_ALLOC(sk, "sk's malloc failed in alloc_glwe_secret_key");

	sk->N = N;
	sk->k = k;

	sk->values = calloc(N * k, sizeof(double));
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
	for (uint64_t j = 0; j < sk->k; j++)
	{
		// TODO: should we forgo the loop and make it a single call to uniform_random_vec?
		CHECK_CALL(uniform_random_vec(nn, glwe_sk_extract_poly(sk, j), 1, nn, nb_bits),
		           "random vector generation failed in key generation");
	}

	return 0;
cleanup:
	return -1;
}

PolyUniv* glwe_sk_extract_poly(GLWESecretKey* sk, uint64_t pos)
{
	assert(pos >= 0 && pos < sk->k);
	return sk->values + sk->N * pos;
}

void delete_glwe_secret_key(GLWESecretKey* sk)
{
	if (!sk) return;
	free(sk->values);
	free(sk);
}

GLWESecretKeyDFT* alloc_glwe_secret_key_dft(uint64_t N, uint64_t k)
{
	uint64_t j           = 0;
	GLWESecretKeyDFT* sk = malloc(sizeof(GLWESecretKeyDFT));
	CHECK_ALLOC(sk, "sk's malloc failed in new_glwe_secret_key.");
	sk->N = N;
	sk->k = k;

	sk->values = malloc(N * k * sizeof(PolyUnivDFT));
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
	return sk_dft->values + sk_dft->N * pos;
}

void delete_glwe_secret_key_dft(GLWESecretKeyDFT* sk_dft)
{
	if (!sk_dft) return;

	free(sk_dft->values);
	free(sk_dft);
}
