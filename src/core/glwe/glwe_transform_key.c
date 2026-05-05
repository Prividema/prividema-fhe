#include "glwe_transform_key.h"

#include <string.h>

#include "glwe_key.h"
#include "univariate_polynomial.h"

void glwe_sk_prepare(const MODULE* module, GLWESecretKeyPrepared* result_dft, const GLWESecretKey* sk)
{
	uint64_t sk_size = sk->nn * sk->k * sizeof(PolyUniv);
	memcpy(result_dft->values_coef, sk->values, sk_size);
	for (uint64_t j = 0; j < sk->k; j++)
		univ_coefs_to_dft(module, glwe_prepared_sk_extract_poly_dft(result_dft, j), glwe_sk_extract_poly(sk, j));
}
