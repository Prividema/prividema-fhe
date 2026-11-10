#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <math.h>
#include <string.h>

#include "bivariate_polynomial.h"
#include "ggsw_arithmetic.h"
#include "ggsw_ciphertext.h"
#include "ggsw_key.h"
#include "ggsw_params.h"
#include "glwe_arithmetic.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "math_utils.h"
#include "univariate_polynomial.h"

// Tests of ggsw_key.h. Parameters: N = 256, k = 2, K = 12, l_a = l_b = 4; GGSW with k_tilde = 2, kappa_tilde = 12,
// 12 rows. Normal noise of stdev 2^-40.

#define CN 256
#define CK 2
#define CKAPPA 12
#define C_LIMBS 12

typedef struct
{
	PvdaBackend* module;
	GLWEParams* params_glwe;
	GGSWParams* params_ggsw;
	GLWESecretKey* sk;
	GLWESecretKeyPrepared* sk_prep;
} Ctx;

static Ctx ctx_new(void)
{
	Ctx c;
	c.module = pvda_new_spqlios_backend(CN);
	PvdaSeed seed;
	for (int i = 0; i < PVDA_SEED_BYTES; ++i) seed.bytes[i] = (uint8_t)(11 * i + 3);
	cr_assert(pvda_backend_set_rng(c.module, PVDA_RNG_CPU, &seed) == 0);
	c.params_glwe = new_glwe_params(CN, CK, CKAPPA, C_LIMBS, ldexp(1.0, -40), NOISE_NORMAL);
	c.params_ggsw = new_ggsw_params(c.params_glwe, CK, CKAPPA, C_LIMBS);
	c.sk          = new_glwe_secret_key(c.params_glwe);
	c.sk_prep     = new_glwe_secret_key_prepared(c.module, c.params_glwe);
	cr_assert(glwe_sk_random_ternary(c.module, c.sk) == 0);
	cr_assert(glwe_sk_prepare(c.module, c.sk_prep, c.sk) == 0);
	return c;
}

static void ctx_delete(Ctx* c)
{
	delete_glwe_secret_key_prepared(c->sk_prep);
	delete_glwe_secret_key(c->sk);
	delete_ggsw_params(c->params_ggsw);
	delete_glwe_params(c->params_glwe);
	pvda_delete_backend(c->module);
}

/** @brief Decrypts ct into TnX. */
static void decrypt_tnx(Ctx* c, const GLWECiphertext* ct, PolyUnivTnX* out)
{
	PolyBiv* dec = new_biv(ct->params);
	cr_assert(glwe_secret_decrypt(c->module, dec, c->sk_prep, ct) == 0);
	cr_assert(biv_to_univ_tnx(ct->params, out, dec, 0) == 0);
	delete_biv(dec);
}

// =============================================
// |                                           |
// |          GLWEGadget -> GGSW KSKs          |
// |                                           |
// =============================================

Test(ggsw_key, ksk_encrypts_minus_each_key_component)
{
	Ctx c                     = ctx_new();
	GGSWCiphertextPrep* ksk[CK] = {0};
	// One pre-allocated, the other allocated by the function
	ksk[1] = new_ggsw_prep(c.module, c.params_ggsw);
	cr_assert(eq(int, glwegadget_to_ggsw_ksk_generate(c.module, ksk, c.params_ggsw, c.sk_prep), 0));

	// GLWE(1/4) times ksk[i] decrypts to -s_i / 4
	GLWECiphertext* ct  = new_glwe(c.params_glwe);
	GLWECiphertext* res = new_glwe(c.params_glwe);
	PolyUnivTnX m[CN] = {UINT64_C(1) << 62}, got[CN];
	cr_assert(glwe_secret_encrypt_tnx(c.module, ct, c.sk_prep, m) == 0);
	for (int i = 0; i < CK; ++i)
	{
		cr_assert_not_null(ksk[i]);
		cr_assert(eq(int, ggsw_external_product(c.module, res, ct, ksk[i]), 0));
		cr_assert(glwe_normalize(c.module, res, res) == 0);
		decrypt_tnx(&c, res, got);
		const PolyUnivZnx* s_i = glwe_sk_extract_poly(c.sk, i);
		for (int p = 0; p < CN; ++p)
			cr_assert(lt(u64, tnx_torus_distance(got[p], (uint64_t)(-s_i[p]) << 62), UINT64_C(1) << 46),
			          "component %d coef %d", i, p);
	}

	delete_glwe(res);
	delete_glwe(ct);
	for (int i = 0; i < CK; ++i) delete_ggsw_prep(ksk[i]);
	ctx_delete(&c);
}

