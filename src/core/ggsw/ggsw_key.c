#include "ggsw_key.h"

#include "glwe_key.h"
#include "logger.h"
#include "rng.h"
#include "spqlios_alias.h"

//! GGSW PART (begin)

PolyUniv** new_ggsw_secret_key_values(uint64_t N, uint64_t k)
{
	PolyUniv** values = malloc(k * sizeof(PolyUniv*));
	if (log_is_null(values, "values's malloc in new_ggsw_secret_key_values") < 0) return NULL;

	for (uint64_t j = 0; j < k; j++) {
		values[j] = calloc(N, sizeof(int64_t));
		if (log_is_null(values[j], "values elements' calloc failed in new_ggsw_secret_key_values") < 0) {
			for (uint64_t t = 0; t < j; t++) free(values[t]);
			free(values);

			return NULL;
		}
	}

	return values;
}

PolyUniv** new_uniform_ggsw_secret_key_values(const MODULE* module, uint64_t k, uint64_t nb_bits)
{
	uint64_t N        = module->nn;
	PolyUniv** values = malloc(k * sizeof(PolyUniv*));
	if (log_is_null(values, "values's malloc failed in new_uniform_ggsw_secret_key_values") < 0) return NULL;

	// Uniform random generation of k Zn[X] polynomials.
	for (uint64_t j = 0; j < k; j++) {
		values[j] = malloc(N * sizeof(double));
		if (log_is_null(values[j], "values elements' calloc failed in new_uniform_ggsw_secret_key_values") < 0) {
			for (uint64_t t = 0; t < j; t++) free(values[t]);
			free(values);

			return NULL;
		}

		if (inplace_uniform_random_vec(N, values[j], 1, N, nb_bits) < 0) {
			log_perror("inplace_uniform_random_vec failed in new_uniform_ggsw_secret_key_values");
			for (uint64_t t = 0; t < j; t++) free(values[t]);
			free(values);

			return NULL;
		}
	}

	return values;
}

void delete_ggsw_secret_key_values(PolyUniv** values, uint64_t k)
{
	for (uint64_t j = 0; j < k; j++) free(values[j]);
	free(values);
}

GGSWSecretKey* new_ggsw_secret_key(uint64_t N, uint64_t k)
{
	GGSWSecretKey* sk = malloc(sizeof(GGSWSecretKey));
	if (log_is_null(sk, "sk's malloc failed in new_ggsw_secret_key") < 0) return NULL;

	sk->N = N;
	sk->k = k;

	sk->values = new_ggsw_secret_key_values(N, k); 
	if (log_is_null(sk->values, "new_ggsw_secret_key_values failed in new_ggsw_secret_key") < 0)
		return NULL;

	return sk;
}

GGSWSecretKey* new_uniform_ggsw_secret_key(const MODULE* module, uint64_t k, uint64_t nb_bits)
{
	uint64_t N        = module->nn;
	GGSWSecretKey* sk = malloc(sizeof(GGSWSecretKey));
	if (log_is_null(sk, "sk's malloc failed in new_uniform_ggsw_secret_key") < 0) return NULL;

	sk->N      = N;
	sk->k      = k;
	sk->values = new_uniform_ggsw_secret_key_values(module, k, nb_bits);

	return sk;
}

void delete_ggsw_secret_key(GGSWSecretKey* sk)
{
	delete_ggsw_secret_key_values(sk->values, sk->k);
	free(sk);
}

GGSWSecretKey* transform_ggsw_secret_key_dft_to_not_dft(const MODULE* module, const GGSWSecretKeyDFT* sk_dft)
{
	GGSWSecretKey* sk = new_ggsw_secret_key(sk_dft->N, sk_dft->k);
	if (sk == NULL)
	{
		log_perror("new_ggsw_secret_key failed in transform_ggsw_secret_key_dft_to_not_dft");
		return NULL;
	}

	for (uint64_t j = 0; j < sk->k; j++)
		if (vec_znx_idft_p(module, sk->values[j], 1, sk_dft->values[j], 1) < 0)
		{
			log_perror("vec_znx_idft_p failed in transform_ggsw_secret_key_dft_to_not_dft");

			for (uint64_t t = 0; t < j; t++) free(sk->values[j]);
			free(sk->values);
			free(sk);

			return NULL;
		}

	return sk;
}

PolyUniv** transform_ggsw_secret_key_values_dft_to_not_dft(const MODULE* module, const PolyUnivDFT** values_dft, uint64_t k)
{
	uint64_t N        = module->nn;

	PolyUniv** values = malloc(k * sizeof(PolyUniv*));
	if (log_is_null(values, "values's calloc failed in transform_ggsw_secret_key_values_dft_to_not_dft") < 0)
		return NULL;

	for (uint64_t j = 0; j < k; j++) {
		values[j] = calloc(N, sizeof(int64_t));
		if (log_is_null(values[j],
		                "values elements' calloc failed in transform_ggsw_secret_key_values_dft_to_not_dft") < 0) {
			for (uint64_t t = 0; t < j; t++) free(values[t]);
			free(values);

			return NULL;
		}
		if (vec_znx_idft_p(module, values[j], 1, values_dft[j], 1) < 0)
		{
			log_perror("vec_znx_idft_p failed in transform_ggsw_secret_key_values_dft_to_not_dft");

			for (uint64_t t = 0; t < j; t++) free(values[t]);
			free(values);

			return NULL;
		}
	}

	return values;
}

GLWESecretKey* transform_ggsw_secret_key_to_glwe_secret_key(const GGSWSecretKey* sk_ggsw)
{
	GLWESecretKey* sk_glwe = new_glwe_secret_key(sk_ggsw->N, sk_ggsw->k);
	if (log_is_null(sk_glwe, "sk_glwe's malloc failed in transform_ggsw_secret_key_to_glwe_secret_key") < 0)
		return NULL;

	for (uint64_t j = 0; j < sk_ggsw->k; j++)
		for (int p = 0; p < sk_ggsw->N; p++) 
			sk_glwe->values[j][p] = sk_ggsw->values[j][p];

	return sk_glwe;
}

//! GGSW PART in the DFT domain (begin)

PolyUnivDFT** new_ggsw_secret_key_values_dft(uint64_t N, uint64_t k)
{
	PolyUnivDFT** values_dft = malloc(k * sizeof(PolyUnivDFT*));
	if (log_is_null(values_dft, "values_dft's malloc failed in new_ggsw_secret_key_values_dft") < 0) return NULL;

	for (uint64_t j = 0; j < k; j++) {
		values_dft[j] = calloc(N, sizeof(double));
		if (log_is_null(values_dft[j], "values_dft elements' calloc failed in new_ggsw_secret_key_values_dft") < 0) {
			for (uint64_t t = 0; t < j; t++) free(values_dft[t]);
			free(values_dft);

			return NULL;
		}
	}

	return values_dft;
}

PolyUnivDFT** new_uniform_ggsw_secret_key_values_dft(const MODULE* module, uint64_t k, uint64_t nb_bits)
{
	uint64_t N               = module->nn;
	PolyUnivDFT** values_dft = malloc(k * sizeof(PolyUnivDFT*));
	if (log_is_null(values_dft, "values_dft's malloc failed in new_uniform_ggsw_secret_key_values_dft") < 0)
		return NULL;

	// Uniform random generation of k Zn[X] polynomials.
	for (uint64_t j = 0; j < k; j++) {
		values_dft[j] = malloc(N * sizeof(double));
		if (log_is_null(values_dft[j], "values_dft elements' calloc failed in new_uniform_ggsw_secret_key_values_dft") <
		    0) {
			for (uint64_t t = 0; t < j; t++) free(values_dft[t]);
			free(values_dft);

			return NULL;
		}

		if (inplace_uniform_random_vec_znx_dft(module, values_dft[j], 1, nb_bits) < 0) {
			log_perror("inplace_uniform_random_vec_znx_dft failed in new_uniform_ggsw_secret_key_values_dft");

			for (uint64_t t = 0; t < j; t++) free(values_dft[t]);
			free(values_dft);

			return NULL;
		}
	}

	return values_dft;
}

void delete_ggsw_secret_key_values_dft(PolyUnivDFT** values, uint64_t k)
{
	for (uint64_t j = 0; j < k; j++) free(values[j]);
	free(values);
}

GGSWSecretKeyDFT* new_ggsw_secret_key_dft(uint64_t N, uint64_t k)
{
	GGSWSecretKeyDFT* sk_dft = malloc(sizeof(GGSWSecretKeyDFT));
	if (log_is_null(sk_dft, "sk_dft's malloc failed in new_ggsw_secret_key_dft") < 0) return NULL;

	sk_dft->N      = N;
	sk_dft->k      = k;
	sk_dft->values = new_ggsw_secret_key_values_dft(N, k);

	return sk_dft;
}

GGSWSecretKeyDFT* new_uniform_ggsw_secret_key_dft(const MODULE* module, uint64_t k, uint64_t nb_bits)
{
	uint64_t N               = module->nn;
	GGSWSecretKeyDFT* sk_dft = malloc(sizeof(GGSWSecretKeyDFT));

	if (log_is_null(sk_dft, "sk_dft's malloc failed in new_uniform_ggsw_secret_key_dft") < 0) return NULL;
	sk_dft->N      = N;
	sk_dft->k      = k;
	sk_dft->values = new_uniform_ggsw_secret_key_values_dft(module, k, nb_bits);

	return sk_dft;
}

void delete_ggsw_secret_key_dft(GGSWSecretKeyDFT* sk_dft)
{
	delete_ggsw_secret_key_values_dft(sk_dft->values, sk_dft->k);
	free(sk_dft);
}

GGSWSecretKeyDFT* transform_ggsw_secret_key_not_dft_to_dft(const MODULE* module, const GGSWSecretKey* sk)
{
	GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(sk->N, sk->k);

	for (uint64_t j = 0; j < sk_dft->k; j++) vec_znx_dft_p(module, sk_dft->values[j], 1, sk->values[j], 1, sk->N);

	return sk_dft;
}

PolyUnivDFT** transform_ggsw_secret_key_values_not_dft_to_dft(const MODULE* module, const PolyUniv** values, uint64_t k)
{
	uint64_t N               = module->nn;

	PolyUnivDFT** values_dft = malloc(k * sizeof(PolyUniv*));
	if (log_is_null(values_dft, "values_dft's malloc failed in transform_ggsw_secret_key_values_not_dft_to_dft") < 0)
		return NULL;

	for (uint64_t j = 0; j < k; j++) {
		values_dft[j] = calloc(N, sizeof(double));
		if (log_is_null(values_dft[j],
		                "values_dft elements' calloc failed in transform_ggsw_secret_key_values_not_dft_to_dft") < 0) {
			for (uint64_t t = 0; t < j; t++) free(values_dft[t]);
			free(values_dft);

			return NULL;
		}

		vec_znx_dft_p(module, values_dft[j], 1, values[j], 1, N);
	}

	return values_dft;
}

GLWESecretKeyDFT* transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(const GGSWSecretKeyDFT* sk_ggsw_dft)
{
	GLWESecretKeyDFT* sk_glwe_dft = new_glwe_secret_key_dft(sk_ggsw_dft->N, sk_ggsw_dft->k);
	if (log_is_null(sk_glwe_dft,
	                "sk_glwe_dft's malloc failed in transform_ggsw_secret_key_dft_to_glwe_secret_key_dft") < 0)
		return NULL;

	for (uint64_t j = 0; j < sk_ggsw_dft->k; j++)
		for (int p = 0; p < sk_ggsw_dft->N; p++) sk_glwe_dft->values[j][p] = sk_ggsw_dft->values[j][p];

	return sk_glwe_dft;
}