#include "ggsw_key.h"

#include "glwe_key.h"
#include "logger.h"
#include "rng.h"
#include "utils.h"
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

int uniform_ggsw_secret_key(const MODULE* module, GGSWSecretKey* sk, uint64_t nb_bits)
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

void delete_ggsw_secret_key(GGSWSecretKey* sk)
{
	delete_ggsw_secret_key_values(sk->values, sk->k);
	free(sk);
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

int uniform_ggsw_secret_key_dft(const MODULE* module, GGSWSecretKeyDFT* sk_dft, uint64_t nb_bits)
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

void delete_ggsw_secret_key_dft(GGSWSecretKeyDFT* sk_dft)
{
	delete_ggsw_secret_key_values_dft(sk_dft->values, sk_dft->k);
	free(sk_dft);
}
