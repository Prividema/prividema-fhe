#include "ggsw_key.h"

#include <stdint.h>
#include <stdlib.h>

#include "glwe_key.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "univariate_polynomial.h"
#include "utils.h"

int generate_glwegad_to_ggsw_ksk(const MODULE* module, GGSWCiphertextPrep** ggsw_ksks, const GGSWParams* ggsw_params,
                                 const GLWESecretKeyPrepared* sk_prep)
{
	int status  = -1;
	uint64_t nn = ggsw_params->params_glwe->nn;
	uint64_t k  = ggsw_params->params_glwe->k;

	PolyUniv* neg_sk_i       = NULL;
	GGSWCiphertext* tmp_ggsw = NULL;
	neg_sk_i                 = new_univ(ggsw_params->params_glwe);
	CHECK_ALLOC(neg_sk_i, "allocation failed in ggsw ksk generation");

	tmp_ggsw = new_ggsw(ggsw_params);
	CHECK_ALLOC(tmp_ggsw, "allocation failed in ggsw ksk generation");
	for (uint64_t i = 0; i < k; ++i)
	{
		if (!ggsw_ksks[i])
		{
			ggsw_ksks[i] = new_ggsw_prep(ggsw_params);
			CHECK_ALLOC(ggsw_ksks[i], "Destination GGSW allocation failed in GGSW KSK generation");
		}
		for (int p = 0; p < nn; ++p)
		{
			neg_sk_i[p] = -glwe_prepared_sk_extract_poly_coefs(sk_prep, i)[p];
		}
		CHECK_CALL(ggsw_secret_encrypt(module, tmp_ggsw, sk_prep, neg_sk_i),
		           "GGSW encryption failed in GGSW KSK generation");
		CHECK_CALL(ggsw_prepare(module, ggsw_ksks[i], tmp_ggsw), "GGSW preparation failed in GGSW KSK generation");
	}

	status = 0;
cleanup:
	delete_ggsw(tmp_ggsw);
	delete_univ(neg_sk_i);
	return status;
}

//TODO: add public key crypto
