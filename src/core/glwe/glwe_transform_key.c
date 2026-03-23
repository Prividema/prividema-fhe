#include "glwe_transform_key.h"

#include "glwe_key.h"

void transform_glwe_secret_key_not_dft_to_dft(const MODULE* module, GLWESecretKeyDFT* result_dft,
                                              const GLWESecretKey* sk)
{
	for (uint64_t j = 0; j < sk->k; j++)
		pvda_vec_znx_dft(module, glwe_sk_dft_retrieve_vec_pos(result_dft, j), 1, glwe_sk_retrieve_vec_pos(sk, j), 1,
		                 sk->N);
}
