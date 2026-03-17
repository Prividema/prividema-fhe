#include "glwe_key.h"

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

	sk->values = calloc(k, sizeof(double*));
	CHECK_ALLOC(sk->values, "values creation failed in glwe key generation");

	for (j = 0; j < k; j++)
	{
		sk->values[j] = calloc(N, sizeof(double));
		CHECK_ALLOC(sk->values[j], "values elements' calloc failed in alloc_glwe_secret_key");
	}

	return sk;
cleanup:
	for (uint64_t t = 0; t < j; t++) free(sk->values[t]);
	if (sk) free(sk->values);
	free(sk);
	return NULL;
}

int uniform_glwe_secret_key(const MODULE* module, GLWESecretKey* sk, uint64_t nb_bits)
{
	uint64_t N = module->nn;
	// The Secret key values
	// Uniform random generation of k Zn[X] polynomials.
	for (uint64_t j = 0; j < sk->k; j++)
	{
		CHECK_CALL(uniform_random_vec(N, sk->values[j], 1, N, nb_bits),
		           "random vector generation failed in key generation");
	}

	return 0;
cleanup:
	return -1;
}

void delete_glwe_secret_key(GLWESecretKey* sk)
{
	if (!sk) return;
	for (uint64_t j = 0; j < sk->k; j++) free(sk->values[j]);
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

	sk->values = malloc(k * sizeof(PolyUnivDFT*));
	CHECK_ALLOC(sk->values, "values' malloc failed in alloc_glwe_secret_key_dft");

	for (j = 0; j < k; j++)
	{
		sk->values[j] = calloc(N, sizeof(double));
		CHECK_ALLOC(sk->values[j], "values elements' calloc failed in alloc_glwe_secret_key_dft");
	}

	return sk;
cleanup:

	for (uint64_t t = 0; t < j; t++) free(sk->values[t]);
	if (sk) free(sk->values);
	free(sk);
	return NULL;
}

void delete_glwe_secret_key_dft(GLWESecretKeyDFT* sk_dft)
{
	if (!sk_dft) return;

	for (uint64_t j = 0; j < sk_dft->k; j++) free(sk_dft->values[j]);
	free(sk_dft->values);
	free(sk_dft);
}
