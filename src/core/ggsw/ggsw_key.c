#include "ggsw_key.h"

#include <stdint.h>
#include <stdlib.h>

#include "glwe_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

int glwegadget_to_ggsw_ksk_generate(const PvdaBackend* module, GGSWCiphertextPrep** ggsw_ksks,
                                    const GGSWParams* params_ggsw, const GLWESecretKeyPrepared* sk_prep)
{
	int status  = -1;
	uint64_t nn = params_ggsw->params_glwe->nn;
	uint64_t k  = params_ggsw->params_glwe->k;

	PolyUniv* neg_sk_i       = NULL;
	GGSWCiphertext* tmp_ggsw = NULL;
	neg_sk_i                 = new_univ(params_ggsw->params_glwe);
	CHECK_ALLOC(neg_sk_i, "allocation failed in ggsw ksk generation");

	tmp_ggsw = new_ggsw(params_ggsw);
	CHECK_ALLOC(tmp_ggsw, "allocation failed in ggsw ksk generation");
	for (uint64_t i = 0; i < k; ++i)
	{
		if (!ggsw_ksks[i])
		{
			ggsw_ksks[i] = new_ggsw_prep(params_ggsw);
			CHECK_ALLOC(ggsw_ksks[i], "Destination GGSW allocation failed in GGSW KSK generation");
		}
		for (int p = 0; p < nn; ++p)
		{
			neg_sk_i[p] = -glwe_sk_prepared_extract_poly_coefs(sk_prep, i)[p];
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
