#include "glwe_key.h"
#include <stdint.h>
#include "logger.h"
#include "utils.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "utils.h"

//! bivGLWE PART (begin)

PolyUniv** alloc_glwe_secret_key_values(uint64_t N, uint64_t k)
{
	PolyUniv** values = malloc(k * sizeof(PolyUniv*));
	if (log_is_null(values, "values' malloc failed in new_glwe_secret_key_values.") < 0) return NULL;

	for (uint64_t j = 0; j < k; j++)
	{
		values[j] = calloc(N, sizeof(int64_t));
		if (log_is_null(values[j], "values elements' calloc failed in new_glwe_secret_key_values.") < 0)
		{
			for (uint64_t t = 0; t < j; t++) free(values[t]);
			free(values);

			return NULL;
		}
	}

	return values;
}

PolyUnivDFT** alloc_glwe_secret_key_values_dft(uint64_t N, uint64_t k)
{
	PolyUnivDFT** values = malloc(k * sizeof(PolyUnivDFT*));
	if (log_is_null(values, "values' malloc failed in new_glwe_secret_key_values.") < 0) return NULL;

	for (uint64_t j = 0; j < k; j++)
	{
		values[j] = calloc(N, sizeof(double));
		if (log_is_null(values[j], "values elements' calloc failed in new_glwe_secret_key_values.") < 0)
		{
			for (uint64_t t = 0; t < j; t++) free(values[t]);
			free(values);

			return NULL;
		}
	}

	return values;
}


void delete_glwe_secret_key_values(PolyUniv** values, uint64_t k)
{
	for (uint64_t j = 0; j < k; j++) free(values[j]);
	free(values);
}

GLWESecretKey* alloc_glwe_secret_key(uint64_t N, uint64_t k)
{
	GLWESecretKey* sk = malloc(sizeof(GLWESecretKey));
	if (log_is_null(sk, "sk's malloc failed in new_glwe_secret_key.") < 0) return NULL;

	sk->N = N;
	sk->k = k;

	sk->values = alloc_glwe_secret_key_values(N, k);
	if (log_is_null(sk->values, "new_glwe_secret_key_values failed in new_glwe_secret_key") < 0)
	{
		free(sk);
		return NULL;
	}

	return sk;
}

int uniform_glwe_secret_key(const MODULE* module, GLWESecretKey* sk, uint64_t nb_bits)
{
	int status = -1;

	// Uniform random generation of k Zn[X] polynomials.
	for (uint64_t j = 0; j < sk->k; j++)
		CHECK_CALL(uniform_random_pol_znx(sk->values[j], sk->N, nb_bits), 
				  "uniform_random_vec failed in new_uniform_ggsw_secret_key_values");		

	status = 0;

cleanup:
	return status;
}

void delete_glwe_secret_key(GLWESecretKey* sk)
{
	delete_glwe_secret_key_values(sk->values, sk->k);
	free(sk);
}


GLWESecretKeyDFT* alloc_glwe_secret_key_dft(uint64_t N, uint64_t k)
{
	GLWESecretKeyDFT* sk = malloc(sizeof(GLWESecretKeyDFT));
	CHECK_ALLOC(sk, "sk's malloc failed in new_glwe_secret_key.");
	sk->N = N;
	sk->k = k;

	sk->values = alloc_glwe_secret_key_values_dft(N, k);

cleanup:
 
  return sk;
}


int uniform_glwe_secret_key_dft(const MODULE* module, GLWESecretKeyDFT* sk_dft, uint64_t nb_bits)
{
	int status = -1;

	// Uniform random generation of k Zn[X] polynomials.
	for (uint64_t j = 0; j < sk_dft->k ; j++)
		CHECK_CALL(uniform_random_vec_znx_dft(module, sk_dft->values[j], 1, nb_bits),
			      "uniform_random_vec_znx_dft failed in new_uniform_ggsw_secret_key_values_dft");

	status = 0;

cleanup:

	return status;
}

void delete_glwe_secret_key_dft(GLWESecretKeyDFT* sk_dft)
{
  if (!sk_dft) return;

	for (uint64_t j = 0; j < sk_dft->k; j++) free(sk_dft->values[j]);
	free(sk_dft->values);
	free(sk_dft);
}
