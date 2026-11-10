#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <criterion/redirect.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "bivariate_polynomial.h"
#include "ggsw_arithmetic.h"
#include "ggsw_ciphertext.h"
#include "ggsw_key.h"
#include "ggsw_params.h"
#include "glwe_arithmetic.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwegadget_arithmetic.h"
#include "glwegadget_ciphertext.h"
#include "glwegadget_key.h"
#include "logger.h"
#include "math_utils.h"
#include "random_polynomial.h"
#include "univariate_polynomial.h"

// Tests of the GGSW / GLWEGadget layer not covered elsewhere: parameters and sizes, views, prepared external
// product, automorphism key collections, single prepared trace expansion and printing.
//
// Parameters: N = 256, k = 2, K = 12, l_a = l_b = 4; GGSW with k_tilde = 2, kappa_tilde = 12, 12 rows
// (l_tilde_a = l_tilde_b = 4). Normal noise of stdev 2^-40.

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
// |            Parameters and sizes           |
// |                                           |
// =============================================

Test(ggsw_params, derived_sizes)
{
	GLWEParams* params_glwe = new_glwe_params(CN, CK, CKAPPA, C_LIMBS, 0, NOISE_UNIFORM_POWER_OF_TWO);

	// 3 * 4 - 1 = 11 rows: l_tilde_a = 4, l_tilde_b = 3
	GGSWParams* p = new_ggsw_params(params_glwe, CK, 7, 11);
	cr_assert_not_null(p);
	cr_assert(eq(ptr, (void*)p->params_glwe, params_glwe));
	cr_assert(eq(u64, p->k_tilde, CK));
	cr_assert(eq(u64, p->kappa_tilde, 7));
	cr_assert(eq(u64, ggsw_num_rows(p), 11));
	cr_assert(eq(u64, ggsw_params_l_tilde_a(p), 4));
	cr_assert(eq(u64, ggsw_params_l_tilde_b(p), 3));
	cr_assert(eq(u64, ggsw_total_n_glwe_limbs(p), 11 * C_LIMBS));
	cr_assert(eq(u64, ggsw_coef_number(p), 11 * C_LIMBS * CN));
	cr_assert(eq(u64, ggsw_coef_number_dft(p), 11 * C_LIMBS * CN / 2));
	cr_assert(eq(u64, ggsw_bytes(p), 11 * C_LIMBS * CN * sizeof(int64_t)));
	delete_ggsw_params(p);
	delete_ggsw_params(NULL);

	GLWEGadgetParams* g = new_glwegadget_params(params_glwe, 9, 5);
	cr_assert_not_null(g);
	cr_assert(eq(ptr, (void*)g->params_glwe, params_glwe));
	cr_assert(eq(u64, g->kappa_tilde, 9));
	cr_assert(eq(u64, g->l_tilde, 5));
	cr_assert(eq(u64, glwegadget_coef_number(g), 5 * glwe_coef_number(params_glwe)));
	delete_glwegadget_params(g);
	delete_glwegadget_params(NULL);

	delete_glwe_params(params_glwe);
}

// =============================================
// |                                           |
// |                   Views                   |
// |                                           |
// =============================================

Test(ggsw_views, glwe_view_and_as_biv)
{
	Ctx c                = ctx_new();
	GGSWCiphertext* ggsw = new_ggsw(c.params_ggsw);
	uint64_t glwe_coefs  = glwe_coef_number(c.params_glwe);

	// Row (prec_lvl, sk_idx) is GLWE number (prec_lvl - 1) * (k_tilde + 1) + sk_idx
	for (int64_t lvl = 1; lvl <= 4; ++lvl)
		for (int64_t j = 0; j <= CK; ++j)
		{
			GLWECiphertext v = ggsw_glwe_view(ggsw, j, lvl);
			cr_assert(eq(ptr, (void*)v.params, c.params_glwe));
			cr_assert(eq(ptr, v.vec, ggsw->mat + ((lvl - 1) * (CK + 1) + j) * glwe_coefs));
		}

	PolyBiv flat = ggsw_as_biv(ggsw);
	cr_assert(eq(u64, flat.nn, CN));
	cr_assert(eq(u64, flat.nb_limbs, ggsw_total_n_glwe_limbs(c.params_ggsw)));
	cr_assert(eq(i64, flat.stride, CN));
	cr_assert(eq(ptr, flat.ptr, ggsw->mat));
	for (uint64_t i = 0; i < ggsw_coef_number(c.params_ggsw); ++i) cr_assert(eq(i64, ggsw->mat[i], 0));

	delete_ggsw(ggsw);
	delete_ggsw(NULL);
	ctx_delete(&c);
}

Test(glwegadget_views, glwe_view)
{
	Ctx c                     = ctx_new();
	GLWEGadgetParams* g       = new_glwegadget_params(c.params_glwe, CKAPPA, 3);
	GLWEGadgetCiphertext* gad = new_glwegadget(g);
	cr_assert_not_null(gad);
	for (uint64_t i = 0; i < glwegadget_coef_number(g); ++i) cr_assert(eq(i64, gad->mat[i], 0));
	for (uint64_t i = 1; i <= 3; ++i)
	{
		GLWECiphertext v = glwegadget_glwe_view(gad, i);
		cr_assert(eq(ptr, (void*)v.params, c.params_glwe));
		cr_assert(eq(ptr, v.vec, gad->mat + (i - 1) * glwe_coef_number(c.params_glwe)));
	}
	delete_glwegadget(gad);
	delete_glwegadget(NULL);
	delete_glwegadget_params(g);
	ctx_delete(&c);
}

// =============================================
// |                                           |
// |             Encryption checks             |
// |                                           |
// =============================================

Test(ggsw_encrypt, rejects_k_tilde_greater_than_k)
{
	set_log_flag(0);
	Ctx c                = ctx_new();
	GGSWParams* bad      = new_ggsw_params(c.params_glwe, CK + 1, CKAPPA, C_LIMBS);
	GGSWCiphertext* ggsw = new_ggsw(bad);
	PolyUnivZnx* m       = new_univ_znx(c.params_glwe);
	memset(m, 0, univ_znx_bytes(c.params_glwe));
	cr_assert(eq(int, ggsw_secret_encrypt(c.module, ggsw, c.sk_prep, m), -1));
	delete_univ_znx(m);
	delete_ggsw(ggsw);
	delete_ggsw_params(bad);
	ctx_delete(&c);
}

// =============================================
// |                                           |
// |         Prepared external product         |
// |                                           |
// =============================================

Test(ggsw_external_product, prepared_matches_unprepared_and_dft_output)
{
	Ctx c                      = ctx_new();
	GGSWCiphertext* ggsw       = new_ggsw(c.params_ggsw);
	GGSWCiphertextPrep* prep   = new_ggsw_prep(c.module, c.params_ggsw);
	GLWECiphertext* ct         = new_glwe(c.params_glwe);
	GLWECiphertext* r_unprep   = new_glwe(c.params_glwe);
	GLWECiphertext* r_prep     = new_glwe(c.params_glwe);
	GLWECiphertext* r_from_dft = new_glwe(c.params_glwe);
	GLWECiphertextDFT* r_dft   = new_glwe_dft(c.module, c.params_glwe);
	PolyUnivZnx* u             = new_univ_znx(c.params_glwe);
	PolyUnivTnX m[CN], got[CN], expected[CN];
	cr_assert_not_null(prep);

	// u = X^5 - 2: u * m is easy to compute exactly on the torus
	memset(u, 0, univ_znx_bytes(c.params_glwe));
	u[0] = -2;
	u[5] = 1;
	cr_assert(univ_tnx_random_uniform(c.module, m, CN) == 0);
	for (int p = 0; p < CN; ++p) expected[p] = -2 * m[p] + (p >= 5 ? m[p - 5] : -m[p - 5 + CN]);

	cr_assert(ggsw_secret_encrypt(c.module, ggsw, c.sk_prep, u) == 0);
	cr_assert(eq(int, ggsw_prepare(c.module, prep, ggsw), 0));
	cr_assert(glwe_secret_encrypt_tnx(c.module, ct, c.sk_prep, m) == 0);

	cr_assert(eq(int, ggsw_external_product_unprepared(c.module, r_unprep, ct, ggsw), 0));
	cr_assert(eq(int, ggsw_external_product(c.module, r_prep, ct, prep), 0));
	cr_assert(eq(int, ggsw_external_product_to_dft(c.module, r_dft, ct, prep), 0));
	cr_assert(eq(int, glwe_dft_to_coef(c.module, r_from_dft, r_dft), 0));

	// The three paths compute the same ciphertext
	cr_assert(memcmp(r_prep->vec, r_unprep->vec, glwe_params_bytes(c.params_glwe)) == 0);
	cr_assert(memcmp(r_prep->vec, r_from_dft->vec, glwe_params_bytes(c.params_glwe)) == 0);

	// And it decrypts to u * m (noise: 12 rows * N coefficients * 2^11 digits * 2^-40 ~ 2^-23)
	cr_assert(glwe_normalize(c.module, r_prep, r_prep) == 0);
	decrypt_tnx(&c, r_prep, got);
	for (int p = 0; p < CN; ++p)
		cr_assert(lt(u64, tnx_torus_distance(got[p], expected[p]), UINT64_C(1) << 46), "coef %d", p);

	delete_univ_znx(u);
	delete_glwe_dft(r_dft);
	delete_glwe(r_from_dft);
	delete_glwe(r_prep);
	delete_glwe(r_unprep);
	delete_glwe(ct);
	delete_ggsw_prep(prep);
	delete_ggsw_prep(NULL);
	delete_ggsw(ggsw);
	ctx_delete(&c);
}

// =============================================
// |                                           |
// |        Automorphism key collections       |
// |                                           |
// =============================================

Test(glwegadget_key, collection_put_get_delete)
{
	Ctx c                                     = ctx_new();
	GLWEGadgetParams* g                       = new_glwegadget_params(c.params_glwe, CKAPPA, 2);
	GLWEAutomorphismKeyCollection* collection = new_automorphism_key_collection(4);
	cr_assert_not_null(collection);
	cr_assert(eq(u64, collection->size, 4));
	for (uint64_t i = 0; i < 4; ++i) cr_assert_null(glwegadget_key_collection_get_key(collection, i));

	GLWEAutomorphismKey* k1 = new_automorphism_key(c.module, g);
	GLWEAutomorphismKey* k2 = new_automorphism_key(c.module, g);
	cr_assert_not_null(k1);
	cr_assert(eq(int, k1->automorphism_p, -2));  // not generated yet
	cr_assert_not_null(k1->glwe_ksk);
	cr_assert(eq(ptr, (void*)k1->glwe_ksk->params, g));

	// Put returns the previous occupant
	cr_assert_null(glwegadget_key_collection_put_key(collection, k1, 1));
	cr_assert(eq(ptr, glwegadget_key_collection_put_key(collection, k2, 1), k1));
	cr_assert(eq(ptr, glwegadget_key_collection_get_key(collection, 1), k2));
	cr_assert_null(glwegadget_key_collection_put_key(collection, k1, 3));
	cr_assert(eq(ptr, glwegadget_key_collection_get_key(collection, 3), k1));

	// Out of range: get gives NULL, put refuses (sentinel -1) and does not store the key
	cr_assert_null(glwegadget_key_collection_get_key(collection, 4));
	cr_assert(eq(ptr, glwegadget_key_collection_put_key(collection, k1, 4), (void*)-1));

	// Deleting with the keys returns how many were deleted
	cr_assert(eq(u64, delete_automorphism_key_collection(collection, 1), 2));

	// Deleting without the keys deletes none of them
	GLWEAutomorphismKeyCollection* other = new_automorphism_key_collection(2);
	GLWEAutomorphismKey* k3              = new_automorphism_key(c.module, g);
	glwegadget_key_collection_put_key(other, k3, 0);
	cr_assert(eq(u64, delete_automorphism_key_collection(other, 0), 0));
	delete_automorphism_key(k3);
	delete_automorphism_key(NULL);

	// Empty collection
	cr_assert(eq(u64, delete_automorphism_key_collection(new_automorphism_key_collection(0), 1), 0));
	cr_assert(eq(u64, delete_automorphism_key_collection(NULL, 1), 0));

	delete_glwegadget_params(g);
	ctx_delete(&c);
}

Test(glwegadget_key, ksk_allocation)
{
	Ctx c               = ctx_new();
	GLWEGadgetParams* g = new_glwegadget_params(c.params_glwe, CKAPPA, 2);
	GLWEKSK* ksk        = new_glwe_ksk(c.module, g);
	cr_assert_not_null(ksk);
	cr_assert(eq(ptr, (void*)ksk->params, g));
	for (int i = 0; i < CK; ++i) cr_assert_not_null(ksk->enc_s[i]);
	delete_glwe_ksk(ksk);
	delete_glwe_ksk(NULL);
	delete_glwegadget_params(g);
	ctx_delete(&c);
}

// =============================================
// |                                           |
// |     Single prepared trace expansion       |
// |                                           |
// =============================================

Test(glwegadget_trace_expand, prepared_single_matches_separate_gadgets)
{
	// k = 1 (trace expansion keys), N = 256, K = 19, l_a = 15; gadget l_tilde = 2
	PvdaBackend* module = pvda_new_spqlios_backend(256);
	PvdaSeed seed       = {{1, 2, 3}};
	cr_assert(pvda_backend_set_rng(module, PVDA_RNG_CPU, &seed) == 0);
	GLWEParams* params_glwe        = new_glwe_params(256, 1, 19, 30, 0, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* g            = new_glwegadget_params(params_glwe, 19, 2);
	GLWESecretKey* sk              = new_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = new_glwe_secret_key_prepared(module, params_glwe);
	cr_assert(glwe_sk_random_binary(module, sk) == 0);
	cr_assert(glwe_sk_prepare(module, sk_prep, sk) == 0);

	GLWEAutomorphismKeyCollection* keys = new_automorphism_key_collection(2 * 256);
	for (uint64_t i = 1; (1ULL << i) <= 256; ++i)
	{
		int64_t p                = 256 / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* key = new_automorphism_key(module, g);
		cr_assert(eq(int, glwe_automorphism_key_generate(module, key, sk_prep, (int)p), 0));
		glwegadget_key_collection_put_key(keys, key, p);
	}

	enum { BUND = 3 };
	PolyUnivZnx* m      = new_univ_znx(params_glwe);
	GLWECiphertext* ct  = new_glwe(params_glwe);
	memset(m, 0, univ_znx_bytes(params_glwe));
	m[0] = 1, m[1] = -1, m[2] = 1;
	cr_assert(glwegadget_packed_secret_encrypt(module, ct, g, sk_prep, m, BUND) == 0);

	// Separate prepared gadgets, and one gadget of depth BUND * l_tilde
	GLWEGadgetCiphertextPrep* sep[BUND] = {0};
	sep[0]                              = new_glwegadget_prep(module, g);
	cr_assert(eq(int, packed_glwegadget_trace_expand_prepared(module, sep, BUND, 2, ct, keys), 0));
	GLWEGadgetParams* big_params    = new_glwegadget_params(params_glwe, 19, BUND * 2);
	GLWEGadgetCiphertextPrep* single = new_glwegadget_prep(module, big_params);
	cr_assert(eq(int, packed_glwegadget_trace_expand_prepared_single(module, single, g, BUND, 2, ct, keys), 0));

	// Half product with a vector that is non-zero only on block r selects gadget r
	PolyBiv* a_small    = new_biv_custom_params(256, 2);
	PolyBiv* a_big      = new_biv_custom_params(256, BUND * 2);
	GLWECiphertext* r1  = new_glwe(params_glwe);
	GLWECiphertext* r2  = new_glwe(params_glwe);
	cr_assert(biv_random_uniform_pow2(module, a_small, 19) == 0);
	for (int r = 0; r < BUND; ++r)
	{
		memset(a_big->ptr, 0, BUND * 2 * 256 * sizeof(int64_t));
		memcpy(a_big->ptr + r * 2 * 256, a_small->ptr, 2 * 256 * sizeof(int64_t));
		cr_assert(eq(int, glwegadget_half_prod(module, r1, sep[r], a_small), 0));
		cr_assert(eq(int, glwegadget_half_prod(module, r2, single, a_big), 0));
		cr_assert(memcmp(r1->vec, r2->vec, glwe_params_bytes(params_glwe)) == 0, "gadget %d", r);
	}

	delete_glwe(r2);
	delete_glwe(r1);
	delete_biv(a_big);
	delete_biv(a_small);
	delete_glwegadget_prep(single);
	delete_glwegadget_params(big_params);
	for (int r = 0; r < BUND; ++r) delete_glwegadget_prep(sep[r]);
	delete_glwe(ct);
	delete_univ_znx(m);
	delete_automorphism_key_collection(keys, 1);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_glwe_secret_key(sk);
	delete_glwegadget_params(g);
	delete_glwe_params(params_glwe);
	pvda_delete_backend(module);
}

// =============================================
// |                                           |
// |                   Print                   |
// |                                           |
// =============================================

Test(glwegadget_print, print_coefs_gad, .init = cr_redirect_stdout)
{
	Ctx c                     = ctx_new();
	GLWEGadgetParams* g       = new_glwegadget_params(c.params_glwe, CKAPPA, 2);
	GLWEGadgetCiphertext* gad = new_glwegadget(g);
	// The zero gadget decrypts to 0 at every level
	cr_assert(eq(int, print_coefs_gad(c.module, gad, c.sk_prep, 2), 0));
	fflush(stdout);
	cr_assert_stdout_eq_str("Gadget lvl 000: 0 (0) 0 (0) \nGadget lvl 001: 0 (0) 0 (0) \n");
	delete_glwegadget(gad);
	delete_glwegadget_params(g);
	ctx_delete(&c);
}
