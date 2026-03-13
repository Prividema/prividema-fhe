#include "ggsw_transform_key.h"

#include "utils.h"

//! bivGLWE KEY PART (begin)

int transform_ggsw_secret_key_dft_to_not_dft(const MODULE* module, GGSWSecretKey* result, const GGSWSecretKeyDFT* sk_dft)
{
	int status = -1;

	for (uint64_t j = 0; j < sk_dft->k; j++)
		CHECK_CALL(vec_znx_idft_p(module, result->values[j], 1, sk_dft->values[j], 1), 
				  "vec_znx_idft_p failed in transform_ggsw_secret_key_dft_to_not_dft");

	status = 0;

cleanup:

	return status;
}


void transform_ggsw_secret_key_to_glwe_secret_key(GLWESecretKey* result, const GGSWSecretKey* sk_ggsw)
{
	for (uint64_t j = 0; j < sk_ggsw->k; j++)
		for (int p = 0; p < sk_ggsw->N; p++) 
			result->values[j][p] = sk_ggsw->values[j][p];
}


//! bivGLWE IN DFT SPACE PART (begin)

void transform_ggsw_secret_key_not_dft_to_dft(const MODULE* module, GGSWSecretKeyDFT* result_dft, const GGSWSecretKey* sk)
{
	for (uint64_t j = 0; j < sk->k; j++) 
		vec_znx_dft_p(module, result_dft->values[j], 1, sk->values[j], 1, sk->N);
}

void transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(GLWESecretKeyDFT* result_dft, const GGSWSecretKeyDFT* sk_ggsw_dft)
{
	for (uint64_t j = 0; j < sk_ggsw_dft->k; j++)
		for (int p = 0; p < sk_ggsw_dft->N; p++) 
			result_dft->values[j][p] = sk_ggsw_dft->values[j][p];
}
