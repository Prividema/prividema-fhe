#include "glwe_transform_key.h"

#include "glwe_key.h"
#include "univariate_polynomial.h"

void transform_glwe_secret_key_not_dft_to_dft(const MODULE* module, GLWESecretKeyDFT* result_dft,
                                              const GLWESecretKey* sk)
{
	for (uint64_t j = 0; j < sk->k; j++)
		univ_coefs_to_dft(module, glwe_sk_extract_poly_dft(result_dft, j), glwe_sk_extract_poly(sk, j));
}
