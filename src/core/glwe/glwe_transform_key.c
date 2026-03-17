#include "glwe_transform_key.h"

void transform_glwe_secret_key_not_dft_to_dft(const MODULE* module, GLWESecretKeyDFT* result_dft,
                                              const GLWESecretKey* sk)
{
	for (uint64_t j = 0; j < sk->k; j++) pvda_vec_znx_dft(module, result_dft->values[j], 1, sk->values[j], 1, sk->N);
}
