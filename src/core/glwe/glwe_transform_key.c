#include "glwe_transform_key.h"

#include "utils.h"

int transform_glwe_secret_key_dft_to_not_dft(const MODULE* module, GLWESecretKey* result, const GLWESecretKeyDFT* sk_dft)
{
	int status = -1;

	for (uint64_t j = 0; j < sk_dft->k; j++)
		CHECK_CALL(vec_znx_idft_p(module, result->values[j], 1, sk_dft->values[j], 1), 
				  "vec_znx_idft_p failed in transform_glwe_secret_key_dft_to_not_dft");

	status = 0;

cleanup:

	return status;
}

void transform_glwe_secret_key_not_dft_to_dft(const MODULE* module, GLWESecretKeyDFT* result_dft, const GLWESecretKey* sk)
{
	for (uint64_t j = 0; j < sk->k; j++) 
		vec_znx_dft_p(module, result_dft->values[j], 1, sk->values[j], 1, sk->N);
}




