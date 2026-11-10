#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <criterion/redirect.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "backend_arithmetic.h"
#include "bivariate_polynomial.h"
#include "glwe_arithmetic.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "logger.h"
#include "math_utils.h"
#include "random_polynomial.h"
#include "univariate_polynomial.h"

// Tests of glwe_ciphertext.h (views, encryption, DFT conversions, print) and glwe_arithmetic.h, checked on
// decrypted values. Parameters: N = 256, k = 2, K = 12, l_a = l_b = 4 (48 bits), normal noise of stdev 2^-40.

#define GN 256
#define GK 2
#define GKAPPA 12
#define G_LIMBS 12

// Decryption error bound: far above the noise (2^-40 per encryption, a few encryptions summed) and far below the
// messages (uniform on the torus): 2^-30, i.e. 2^34 in TnX units.
#define MAX_ERR (UINT64_C(1) << 34)

typedef struct
{
	PvdaBackend* module;
	GLWEParams* params;
	GLWESecretKey* sk;
	GLWESecretKeyPrepared* sk_prep;
} Ctx;

static Ctx ctx_new(void)
{
	Ctx c;
	c.module = pvda_new_spqlios_backend(GN);
	PvdaSeed seed;
	for (int i = 0; i < PVDA_SEED_BYTES; ++i) seed.bytes[i] = (uint8_t)(5 * i + 1);
	cr_assert(pvda_backend_set_rng(c.module, PVDA_RNG_CPU, &seed) == 0);
	c.params  = new_glwe_params(GN, GK, GKAPPA, G_LIMBS, ldexp(1.0, -40), NOISE_NORMAL);
	c.sk      = new_glwe_secret_key(c.params);
	c.sk_prep = new_glwe_secret_key_prepared(c.module, c.params);
	cr_assert(glwe_sk_random_ternary(c.module, c.sk) == 0);
	cr_assert(glwe_sk_prepare(c.module, c.sk_prep, c.sk) == 0);
	return c;
}

static void ctx_delete(Ctx* c)
{
	delete_glwe_secret_key_prepared(c->sk_prep);
	delete_glwe_secret_key(c->sk);
	delete_glwe_params(c->params);
	pvda_delete_backend(c->module);
}

/** @brief Encrypts a random TnX message, returned in m. */
static GLWECiphertext* encrypt_random(Ctx* c, PolyUnivTnX* m)
{
	GLWECiphertext* ct = new_glwe(c->params);
	cr_assert(univ_tnx_random_uniform(c->module, m, GN) == 0);
	cr_assert(glwe_secret_encrypt_tnx(c->module, ct, c->sk_prep, m) == 0);
	return ct;
}

/** @brief Asserts that ct decrypts to expected (up to the noise). */
static void assert_decrypts_to(Ctx* c, const GLWECiphertext* ct, const PolyUnivTnX* expected)
{
	PolyBiv* dec = new_biv(c->params);
	PolyUnivTnX got[GN];
	cr_assert(glwe_secret_decrypt(c->module, dec, c->sk_prep, ct) == 0);
	cr_assert(biv_to_univ_tnx(c->params, got, dec, 0) == 0);
	for (int p = 0; p < GN; ++p)
		cr_assert(lt(u64, tnx_torus_distance(got[p], expected[p]), MAX_ERR), "coef %d: %lx vs %lx", p, got[p],
		          expected[p]);
	delete_biv(dec);
}

// =============================================
// |                                           |
// |                  Views                    |
// |                                           |
// =============================================

Test(glwe_views, extract_poly_view_and_as_biv)
{
	GLWEParams* params = new_glwe_params(GN, GK, GKAPPA, G_LIMBS - 1, 0, NOISE_UNIFORM_POWER_OF_TWO);  // l_b = l_a - 1
	GLWECiphertext* ct = new_glwe(params);
	uint64_t l_a = glwe_params_l_a(params), l_b = glwe_params_l_b(params);
	cr_assert(eq(u64, l_b + 1, l_a));

	// Limb-major layout: limb i of element j starts at (i * (k + 1) + j) * N
	for (uint64_t j = 0; j <= GK; ++j)
	{
		PolyBiv v = glwe_extract_poly_view(ct, j);
		cr_assert(eq(u64, v.nn, GN));
		cr_assert(eq(u64, v.nb_limbs, j == GK ? l_b : l_a));
		cr_assert(eq(i64, v.stride, (GK + 1) * GN));
		cr_assert(eq(ptr, v.ptr, ct->vec + j * GN));
		v.ptr[(v.nb_limbs - 1) * v.stride + 3] = (int64_t)j + 1;
		cr_assert(eq(i64, ct->vec[((v.nb_limbs - 1) * (GK + 1) + j) * GN + 3], (int64_t)j + 1));
	}

	PolyBiv flat = glwe_as_biv(ct);
	cr_assert(eq(u64, flat.nn, GN));
	cr_assert(eq(u64, flat.nb_limbs, G_LIMBS - 1));
	cr_assert(eq(i64, flat.stride, GN));
	cr_assert(eq(ptr, flat.ptr, ct->vec));

	delete_glwe(ct);
	delete_glwe(NULL);
	delete_glwe_params(params);
}

Test(glwe_views, new_glwe_is_zero)
{
	GLWEParams* params = new_glwe_params(GN, GK, GKAPPA, G_LIMBS, 0, NOISE_UNIFORM_POWER_OF_TWO);
	GLWECiphertext* ct = new_glwe(params);
	cr_assert(eq(ptr, (void*)ct->params, params));
	for (uint64_t i = 0; i < glwe_coef_number(params); ++i) cr_assert(eq(i64, ct->vec[i], 0));
	delete_glwe(ct);
	delete_glwe_params(params);
}

// =============================================
// |                                           |
// |          Encryption / decryption          |
// |                                           |
// =============================================

Test(glwe_encrypt, tnx_and_rnx_messages)
{
	Ctx c = ctx_new();
	PolyUnivTnX m[GN];
	GLWECiphertext* ct = encrypt_random(&c, m);
	assert_decrypts_to(&c, ct, m);

	// RnX message: same values read as reals in [-1/2, 1/2)
	PolyUnivRnX r[GN];
	cr_assert(univ_tnx_to_rnx(c.params, r, m) == 0);
	cr_assert(glwe_secret_encrypt_rnx(c.module, ct, c.sk_prep, r) == 0);
	assert_decrypts_to(&c, ct, m);

	// Two encryptions of the same message differ (fresh mask and noise)
	GLWECiphertext* ct2 = new_glwe(c.params);
	cr_assert(glwe_secret_encrypt_tnx(c.module, ct2, c.sk_prep, m) == 0);
	cr_assert(memcmp(ct->vec, ct2->vec, glwe_params_bytes(c.params)) != 0);

	delete_glwe(ct2);
	delete_glwe(ct);
	ctx_delete(&c);
}

Test(glwe_encrypt, zero_ciphertext_decrypts_to_zero)
{
	Ctx c                 = ctx_new();
	GLWECiphertext* ct    = new_glwe(c.params);
	PolyUnivTnX zero[GN] = {0};
	assert_decrypts_to(&c, ct, zero);
	delete_glwe(ct);
	ctx_delete(&c);
}

// =============================================
// |                                           |
// |                Arithmetic                 |
// |                                           |
// =============================================

Test(glwe_arith, add_sub_negate)
{
	Ctx c = ctx_new();
	PolyUnivTnX m1[GN], m2[GN], expected[GN];
	GLWECiphertext* c1  = encrypt_random(&c, m1);
	GLWECiphertext* c2  = encrypt_random(&c, m2);
	GLWECiphertext* res = new_glwe(c.params);

	cr_assert(eq(int, glwe_add(c.module, res, c1, c2), 0));
	for (int p = 0; p < GN; ++p) expected[p] = m1[p] + m2[p];
	assert_decrypts_to(&c, res, expected);

	cr_assert(eq(int, glwe_sub(c.module, res, c1, c2), 0));
	for (int p = 0; p < GN; ++p) expected[p] = m1[p] - m2[p];
	assert_decrypts_to(&c, res, expected);

	cr_assert(eq(int, glwe_negate(c.module, res, c1), 0));
	for (int p = 0; p < GN; ++p) expected[p] = -m1[p];
	assert_decrypts_to(&c, res, expected);
	// Coefficient-wise negation
	for (uint64_t i = 0; i < glwe_coef_number(c.params); ++i) cr_assert(eq(i64, res->vec[i], -c1->vec[i]));

	// In place: c1 - c1 = 0 exactly
	cr_assert(eq(int, glwe_sub(c.module, c1, c1, c1), 0));
	for (uint64_t i = 0; i < glwe_coef_number(c.params); ++i) cr_assert(eq(i64, c1->vec[i], 0));

	delete_glwe(res);
	delete_glwe(c2);
	delete_glwe(c1);
	ctx_delete(&c);
}

Test(glwe_arith, normalize_after_many_additions)
{
	Ctx c = ctx_new();
	PolyUnivTnX m[GN], expected[GN];
	GLWECiphertext* ct  = encrypt_random(&c, m);
	GLWECiphertext* acc = new_glwe(c.params);

	// 64 additions without normalization: digits grow beyond [-2^(K-1), 2^(K-1))
	for (int r = 0; r < 64; ++r) cr_assert(glwe_add(c.module, acc, acc, ct) == 0);
	for (int p = 0; p < GN; ++p) expected[p] = 64 * m[p];

	cr_assert(eq(int, glwe_normalize(c.module, acc, acc), 0));
	for (uint64_t i = 0; i < glwe_coef_number(c.params); ++i)
		cr_assert(acc->vec[i] >= -(1 << (GKAPPA - 1)) && acc->vec[i] < (1 << (GKAPPA - 1)), "coef %lu", i);
	// 64 noises of stdev 2^-40: still far below the bound
	assert_decrypts_to(&c, acc, expected);

	delete_glwe(acc);
	delete_glwe(ct);
	ctx_delete(&c);
}

Test(glwe_arith, mult_const_by_monomial_and_dft_paths)
{
	Ctx c = ctx_new();
	PolyUnivTnX m[GN], expected[GN];
	GLWECiphertext* ct  = encrypt_random(&c, m);
	GLWECiphertext* res = new_glwe(c.params);
	GLWECiphertext* res2 = new_glwe(c.params);

	// u = -X^3: m * u is a negacyclic rotation with a sign change
	PolyUnivZnx* u          = new_univ_znx(c.params);
	BE_PolyUnivSvpPrep* u_p = new_univ_svp_prep(c.module);
	memset(u, 0, univ_znx_bytes(c.params));
	u[3] = -1;
	cr_assert(univ_znx_to_svp_prep(c.module, u_p, u) == 0);
	for (int p = 0; p < GN; ++p) expected[p] = p >= 3 ? -m[p - 3] : m[p - 3 + GN];

	cr_assert(eq(int, glwe_mult_const(c.module, res, u_p, ct), 0));
	assert_decrypts_to(&c, res, expected);

	// DFT round trip is exact
	GLWECiphertextDFT* d1 = new_glwe_dft(c.module, c.params);
	GLWECiphertextDFT* d2 = new_glwe_dft(c.module, c.params);
	cr_assert_not_null(d1);
	cr_assert(eq(int, glwe_coef_to_dft(c.module, d1, ct), 0));
	cr_assert(eq(int, glwe_dft_to_coef(c.module, res2, d1), 0));
	cr_assert(memcmp(res2->vec, ct->vec, glwe_params_bytes(c.params)) == 0);

	// glwe_dft_mult_const gives the same ciphertext as glwe_mult_const
	cr_assert(eq(int, glwe_dft_mult_const(c.module, d2, u_p, d1), 0));
	cr_assert(eq(int, glwe_dft_to_coef(c.module, res2, d2), 0));
	cr_assert(memcmp(res2->vec, res->vec, glwe_params_bytes(c.params)) == 0);

	// glwe_dft_add gives the same ciphertext as glwe_add
	cr_assert(eq(int, glwe_dft_add(c.module, d2, d1, d1), 0));
	cr_assert(eq(int, glwe_dft_to_coef(c.module, res2, d2), 0));
	cr_assert(eq(int, glwe_add(c.module, res, ct, ct), 0));
	cr_assert(memcmp(res2->vec, res->vec, glwe_params_bytes(c.params)) == 0);

	delete_glwe_dft(d2);
	delete_glwe_dft(d1);
	delete_glwe_dft(NULL);
	delete_univ_svp_prep(c.module, u_p);
	delete_univ_znx(u);
	delete_glwe(res2);
	delete_glwe(res);
	delete_glwe(ct);
	ctx_delete(&c);
}

// =============================================
// |                                           |
// |                   Print                   |
// |                                           |
// =============================================

Test(glwe_print, print_coefs_of_a_known_message, .init = cr_redirect_stdout)
{
	Ctx c = ctx_new();
	// Noiseless encryption of m = (1/2, 1/4, 0, ...)
	GLWEParams* params = new_glwe_params(GN, GK, GKAPPA, G_LIMBS, 0, NOISE_UNIFORM_POWER_OF_TWO);
	GLWECiphertext* ct = new_glwe(params);
	PolyUnivTnX m[GN]  = {UINT64_C(1) << 63, UINT64_C(1) << 62};
	set_log_flag(LOG_ERROR);  // no warning about the zero noise in the captured output
	cr_assert(glwe_secret_encrypt_tnx(c.module, ct, c.sk_prep, m) == 0);

	cr_assert(eq(int, glwe_print_coefs(c.module, ct, c.sk_prep, 3, 0), 0));
	printf("|");
	// shft = 62: the values rounded to 2 bits
	cr_assert(eq(int, glwe_print_coefs(c.module, ct, c.sk_prep, 2, 62), 0));
	fflush(stdout);
	cr_assert_stdout_eq_str("8000000000000000 (-9223372036854775808) 4000000000000000 (4611686018427387904) 0 (0) |"
	                        "8000000000000000 (2)  4000000000000000 (1)  ");

	delete_glwe(ct);
	delete_glwe_params(params);
	ctx_delete(&c);
}
