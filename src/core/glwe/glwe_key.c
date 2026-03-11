#include "glwe_key.h"
#include "logger.h"
#include "rng.h"
#include "spqlios_alias.h"

//! GLWE PART (begin)

PolyUniv** new_glwe_secret_key_values(uint64_t N, uint64_t k)
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

PolyUniv** new_uniform_glwe_secret_key_values(const MODULE* module, uint64_t k, uint64_t nb_bits)
{
	// The
	uint64_t N = module->nn;

	// The Secret key values
	PolyUniv** values = malloc(k * sizeof(PolyUniv*));
	if (log_is_null(values, "values' malloc failed in new_uniform_glwe_secret_key_values.") < 0) return NULL;

	// Uniform random generation of k Zn[X] polynomials.
	for (uint64_t j = 0; j < k; j++)
	{
		values[j] = calloc(N, sizeof(int64_t));
		if (log_is_null(values[j], "values elements calloc failed in new_uniform_glwe_secret_key_values.") < 0)
		{
			for (uint64_t t = 0; t < j; t++) free(values[t]);
			free(values);

			return NULL;
		}

		if (uniform_random_vec(N, values[j], 1, N, nb_bits) < 0)
		{
			log_perror("uniform_random_vec failed in new_uniform_glwe_secret_key_values.");

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

GLWESecretKey* new_glwe_secret_key(uint64_t N, uint64_t k)
{
	GLWESecretKey* sk = malloc(sizeof(GLWESecretKey));
	if (log_is_null(sk, "sk's malloc failed in new_glwe_secret_key.") < 0) return NULL;

	sk->N = N;
	sk->k = k;

	sk->values = new_glwe_secret_key_values(N, k);
	if (log_is_null(sk->values, "new_glwe_secret_key_values failed in new_glwe_secret_key") < 0)
	{
		free(sk);
		return NULL;
	}

	return sk;
}

GLWESecretKey* new_uniform_glwe_secret_key(const MODULE* module, uint64_t k, uint64_t nb_bits)
{
	uint64_t N = module->nn;

	GLWESecretKey* sk = malloc(sizeof(GLWESecretKey));
	if (log_is_null(sk, "sk's malloc failed in new_uniform_glwe_secret_key.") < 0) return NULL;

	sk->N = N;
	sk->k = k;

	sk->values = new_uniform_glwe_secret_key_values(module, k, nb_bits);
	if (log_is_null(sk->values, "new_uniform_glwe_secret_key_values in new_uniform_glwe_secret_key") < 0)
	{
		free(sk);
		return NULL;
	}

	return sk;
}

void delete_glwe_secret_key(GLWESecretKey* sk)
{
	delete_glwe_secret_key_values(sk->values, sk->k);
	free(sk);
}

GLWESecretKey* transform_glwe_secret_key_dft_to_not_dft(const MODULE* module, const GLWESecretKeyDFT* sk_dft)
{
	GLWESecretKey* sk = new_glwe_secret_key(sk_dft->N, sk_dft->k);
	if (log_is_null(sk, "new_glwe_secret_key failed in transform_glwe_secret_key_dft_to_not_dft") < 0) return NULL;

	for (uint64_t j = 0; j < sk->k; j++)
	{
		if (vec_znx_idft_p(module, sk->values[j], 1, sk_dft->values[j], 1) < 0)
		{
			log_perror("vec_znx_idft_p failed in transform_glwe_secret_key_dft_to_not_dft");

			for (uint64_t t = 0 ; t < j ; t++) free(sk->values[j]);
			free(sk->values);
			free(sk);

			return NULL;
		}
	}

	return sk;
}

PolyUniv** transform_glwe_secret_key_values_dft_to_not_dft(const MODULE* module, const PolyUnivDFT** values_dft, uint64_t k)
{
	uint64_t N = module->nn;

	PolyUniv** values = malloc(k * sizeof(PolyUniv*));
	if (log_is_null(values, "values' malloc failed in transform_glwe_secret_key_values_dft_to_not_dft.") < 0)
		return NULL;

	for (uint64_t j = 0; j < k; j++)
	{
		values[j] = calloc(N, sizeof(int64_t));
		if (log_is_null(values[j],
		                "values elements' malloc failed in transform_glwe_secret_key_values_dft_to_not_dft.") < 0)
		{
			for (uint64_t t = 0; t < j; t++) free(values[t]);
			free(values);

			return NULL;
		}
		if (vec_znx_idft_p(module, values[j], 1, values_dft[j], 1) < 0)
		{
			log_perror("vec_znx_idft_p failed in transform_glwe_secret_key_values_dft_to_not_dft");

			for (uint64_t t = 0; t < j; t++) free(values[t]);
			free(values);

			return NULL;
		}
	}

	return values;
}

//! GLWE PART in the DFT domain (begin)

PolyUnivDFT** new_glwe_secret_key_values_dft(uint64_t N, uint64_t k)
{
	PolyUnivDFT** values_dft = malloc(k * sizeof(PolyUnivDFT*));
	if (log_is_null(values_dft, "values malloc failed in new_glwe_secret_key_values_dft.") < 0) return NULL;

	for (uint64_t j = 0; j < k; j++)
	{
		values_dft[j] = calloc(N, sizeof(double));
		if (log_is_null(values_dft[j], "values elements' calloc failed in new_glwe_secret_key_values_dft.") < 0)
		{
			for (uint64_t t = 0; t < j; t++) free(values_dft[t]);
			free(values_dft);

			return NULL;
		}
	}

	return values_dft;
}

PolyUnivDFT** new_uniform_glwe_secret_key_values_dft(const MODULE* module, uint64_t k, uint64_t nb_bits)
{
	uint64_t N = module->nn;

	PolyUnivDFT** values_dft = malloc(k * sizeof(PolyUnivDFT*));
	if (log_is_null(values_dft, "values_dft's malloc failed in new_uniform_glwe_secret_key_values_dft.") < 0)
		return NULL;

	// Uniform random generation of k Zn[X] polynomials.
	for (uint64_t j = 0; j < k; j++)
	{
		values_dft[j] = malloc(N * sizeof(double));
		if (log_is_null(values_dft[j], "values elements' malloc failed in new_uniform_glwe_secret_key_values_dft.") < 0)
		{
			for (uint64_t t = 0; t < j; t++) free(values_dft[t]);
			free(values_dft);

			return NULL;
		}

		if (uniform_random_vec_znx_dft(module, values_dft[j], 1, nb_bits) < 0)
		{
			log_perror("uniform_random_vec_znx_dft failed in new_uniform_glwe_secret_key_values_dft.");

			for (uint64_t t = 0; t < j; t++) free(values_dft[t]);
			free(values_dft);

			return NULL;
		}
	}

	return values_dft;
}

void delete_glwe_secret_key_values_dft(PolyUnivDFT** values, uint64_t k)
{
	for (uint64_t j = 0; j < k; j++) free(values[j]);

	free(values);
}

GLWESecretKeyDFT* new_glwe_secret_key_dft(uint64_t N, uint64_t k)
{
	GLWESecretKeyDFT* sk_dft = malloc(sizeof(GLWESecretKeyDFT));
	if (log_is_null(sk_dft, "sk_dft's malloc failed in new_glwe_secret_key_dft.") < 0) return NULL;

	sk_dft->N = N;
	sk_dft->k = k;

	sk_dft->values = new_glwe_secret_key_values_dft(N, k);
	if (log_is_null(sk_dft->values, "new_glwe_secret_key_values_dft failed in new_glwe_secret_key_dft") < 0)
	{
		free(sk_dft);
		return NULL;
	}

	return sk_dft;
}

GLWESecretKeyDFT* new_uniform_glwe_secret_key_dft(const MODULE* module, uint64_t k, uint64_t nb_bits)
{
	uint64_t N = module->nn;

	GLWESecretKeyDFT* sk_dft = malloc(sizeof(GLWESecretKeyDFT));
	if (log_is_null(sk_dft, "sk_dft's malloc failed in new_uniform_glwe_secret_key_dft") < 0) return NULL;

	sk_dft->N = N;
	sk_dft->k = k;

	sk_dft->values = new_uniform_glwe_secret_key_values_dft(module, k, nb_bits);
	if (log_is_null(sk_dft->values, "new_uniform_glwe_secret_key_values_dft in new_uniform_glwe_secret_key_dft") < 0)
	{
		free(sk_dft);
		return NULL;
	}
	return sk_dft;
}

void delete_glwe_secret_key_dft(GLWESecretKeyDFT* sk_dft)
{
	delete_glwe_secret_key_values_dft(sk_dft->values, sk_dft->k);
	free(sk_dft);
}

GLWESecretKeyDFT* transform_glwe_secret_key_not_dft_to_dft(const MODULE* module, const GLWESecretKey* sk)
{
	GLWESecretKeyDFT* sk_dft = new_glwe_secret_key_dft(sk->N, sk->k);
	if (log_is_null(sk_dft, "new_glwe_secret_key_dft failed in transform_glwe_secret_key_not_dft_to_dft") < 0)
		return NULL;

	for (uint64_t j = 0; j < sk_dft->k; j++) vec_znx_dft_p(module, sk_dft->values[j], 1, sk->values[j], 1, sk->N);

	return sk_dft;
}

PolyUnivDFT** transform_glwe_secret_key_values_not_dft_to_dft(const MODULE* module, const PolyUniv** values, uint64_t k)
{
	uint64_t N = module->nn;

	PolyUnivDFT** values_dft = malloc(k * sizeof(PolyUnivDFT*));
	if (log_is_null(values_dft, "values_dft malloc failed in transform_glwe_secret_key_values_not_dft_to_dft.") < 0)
		return NULL;

	for (uint64_t j = 0; j < k; j++)
	{
		values_dft[j] = calloc(N, sizeof(double));
		if (log_is_null(values_dft[j],
		                "values_dft elements' calloc failed in transform_glwe_secret_key_values_not_dft_to_dft") < 0)
		{
			for (uint64_t t = 0; t < j; t++) free(values_dft[t]);
			free(values_dft);

			return NULL;
		}
		vec_znx_dft_p(module, values_dft[j], 1, values[j], 1, N);
	}

	return values_dft;
}
