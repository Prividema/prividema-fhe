#include "glwe_transform_key.h"

#include <string.h>

#include "glwe_key.h"
#include "univariate_polynomial.h"

void glwe_sk_prepare(const MODULE* module, GLWESecretKeyPrepared* result_dft, const GLWESecretKey* sk)
{
	uint64_t sk_size = sk->nn * sk->k * sizeof(PolyUniv);
	memcpy(result_dft->values_coef, sk->values, sk_size);
	pvda_vec_znx_dft(module, result_dft->values, sk->k, sk->values, sk->k, sk->nn);
}
