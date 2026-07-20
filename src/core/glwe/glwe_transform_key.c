#include "glwe_transform_key.h"

#include <string.h>

#include "backend_arithmetic.h"
#include "bivariate_polynomial.h"
#include "glwe_key.h"
#include "maths_structures.h"

void glwe_sk_prepare(const PvdaBackend* module, GLWESecretKeyPrepared* result_dft, const GLWESecretKey* sk)
{
	uint64_t sk_size = sk->nn * sk->k * sizeof(PolyUniv);
	memcpy(result_dft->values_coef, sk->values, sk_size);
	uint64_t nn    = sk->nn;
	PolyBiv sk_biv = new_biv_view(nn, sk->k, nn, sk->values);
	pvda_vec_znx_dft(module, result_dft->values, sk->k, &sk_biv);
}
