#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <stdint.h>
#include <sys/types.h>

#include "backend_private.h"
#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "heint.h"
#include "pvda_ffts.h"
#include "random_polynomial.h"
#include "univariate_polynomial.h"

Test(heint_enc, sanity1)
{
	uint64_t nn             = 8;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 4, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 17;  // 2*nn + 1

	PolyUnivTnX* in_tnx = new_univ_tnx(params_glwe);
	PolyBiv* biv_out    = new_biv(params_glwe);

	ntt_table_generate(module, t);
	univ_znx_random_uniform(module, in_tnx, nn, 0, t - 1);

	int st = heint_encode(module, params_glwe, biv_out, nn, t, in_tnx);

	cr_assert(eq(i32, st, 0));
	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_univ_tnx(in_tnx);
	delete_biv(biv_out);
}

Test(heint_enc, backforth)
{
	uint64_t nn             = 8;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 6, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 17;  // 2*nn + 1

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);
	PolyBiv* biv_out     = new_biv(params_glwe);

	ntt_table_generate(module, t);
	univ_znx_random_uniform(module, in_tnx, nn, 0, t - 1);

	int st  = heint_encode(module, params_glwe, biv_out, nn, t, in_tnx);
	int st2 = heint_decode(module, params_glwe, out_tnx, nn, t, biv_out);

	cr_assert(eq(i32, st, 0));
	cr_assert(eq(i32, st2, 0));

	for (int i = 0; i < nn; ++i)
	{
		cr_assert(eq(i64, in_tnx[i], out_tnx[i]));
	}

	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_biv(biv_out);
}
Test(heint_enc, backforth3)
{
	uint64_t nn             = 1024;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 6, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 4085761;

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);
	PolyBiv* biv_out     = new_biv(params_glwe);

	ntt_table_generate(module, t);
	univ_znx_random_uniform(module, in_tnx, nn, 0, t - 1);

	int st  = heint_encode(module, params_glwe, biv_out, nn, t, in_tnx);
	int st2 = heint_decode(module, params_glwe, out_tnx, nn, t, biv_out);

	cr_assert(eq(i32, st, 0));
	cr_assert(eq(i32, st2, 0));

	for (int i = 0; i < nn; ++i)
	{
		cr_assert(eq(i64, in_tnx[i], out_tnx[i]));
	}

	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_biv(biv_out);
}

Test(heint_enc, backforth2)
{
	uint64_t nn             = 1024;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 6, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 12289;  // 6*2*nn + 1

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);
	PolyBiv* biv_out     = new_biv(params_glwe);

	ntt_table_generate(module, t);
	univ_znx_random_uniform(module, in_tnx, nn, 0, t - 1);

	int st  = heint_encode(module, params_glwe, biv_out, nn, t, in_tnx);
	int st2 = heint_decode(module, params_glwe, out_tnx, nn, t, biv_out);

	cr_assert(eq(i32, st, 0));
	cr_assert(eq(i32, st2, 0));

	for (int i = 0; i < nn; ++i)
	{
		cr_assert(eq(i64, in_tnx[i], out_tnx[i]));
	}

	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_biv(biv_out);
}

// heint.h allows any prime t <= 2^31 with t = 1 mod 2N. With N = 1024 and t close to 2^31, the NTT sums
// (N terms in [0, t)) exceed 2^32, so this checks that they are reduced mod t and not mod 2^32.
#define HEINT_LARGE_T 2147473409ull  // Largest prime < 2^31 with t = 1 mod 2048

Test(heint_enc, backforth_large_t)
{
	uint64_t nn             = 1024;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 6, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = HEINT_LARGE_T;

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);
	PolyBiv* biv_out     = new_biv(params_glwe);

	cr_assert(eq(i32, ntt_table_generate(module, t), 0));
	univ_znx_random_uniform(module, in_tnx, nn, 0, t - 1);

	int st  = heint_encode(module, params_glwe, biv_out, nn, t, in_tnx);
	int st2 = heint_decode(module, params_glwe, out_tnx, nn, t, biv_out);

	cr_assert(eq(i32, st, 0));
	cr_assert(eq(i32, st2, 0));

	for (int i = 0; i < nn; ++i)
	{
		cr_assert(eq(u64, in_tnx[i], out_tnx[i]), "slot %d", i);
	}

	pvda_delete_backend(module);
	delete_glwe_params(params_glwe);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_biv(biv_out);
}

Test(heint_internals, backforth_large_t)
{
	uint64_t nn             = 1024;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 6, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = HEINT_LARGE_T;

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* mid_tnx = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);

	cr_assert(eq(i32, ntt_table_generate(module, t), 0));
	univ_znx_random_uniform(module, in_tnx, nn, 0, t - 1);

	int st  = internal_slow_intt_heint(module, mid_tnx, in_tnx, t);
	int st2 = internal_slow_ntt_heint(module, out_tnx, mid_tnx, t);
	cr_assert(eq(i32, st, 0));
	cr_assert(eq(i32, st2, 0));

	for (int i = 0; i < nn; ++i)
	{
		cr_assert(lt(u64, mid_tnx[i], t), "intt output %d not reduced mod t", i);
		cr_assert(eq(u64, in_tnx[i], out_tnx[i]), "slot %d", i);
	}

	pvda_delete_backend(module);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_univ_tnx(mid_tnx);
	delete_glwe_params(params_glwe);
}

Test(heint_internals, sanity1)
{
	uint64_t nn             = 8;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 4, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 17;  // 2*nn + 1

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);

	ntt_table_generate(module, t);
	univ_znx_random_uniform(module, in_tnx, nn, 0, t - 1);

	int st = internal_slow_intt_heint(module, out_tnx, in_tnx, t);

	cr_assert(eq(i32, st, 0));
	pvda_delete_backend(module);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_glwe_params(params_glwe);
}

Test(heint_internals, sanity2)
{
	uint64_t nn             = 8;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 4, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 17;  // 2*nn + 1

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);

	ntt_table_generate(module, t);
	univ_znx_random_uniform(module, in_tnx, nn, 0, t - 1);

	int st = internal_slow_ntt_heint(module, out_tnx, in_tnx, t);

	cr_assert(eq(i32, st, 0));
	pvda_delete_backend(module);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_glwe_params(params_glwe);
}

Test(heint_internals, backforth)
{
	uint64_t nn             = 8;
	PvdaBackend* module     = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 4, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 17;  // 2*nn + 1

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* mid_tnx = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);

	ntt_table_generate(module, t);
	univ_znx_random_uniform(module, in_tnx, nn, 0, t - 1);

	int st  = internal_slow_intt_heint(module, mid_tnx, in_tnx, t);
	int st2 = internal_slow_ntt_heint(module, out_tnx, mid_tnx, t);
	cr_assert(eq(i32, st, 0));
	cr_assert(eq(i32, st2, 0));

	for (int i = 0; i < nn; ++i)
	{
		cr_assert(eq(i64, in_tnx[i], out_tnx[i]));
	}

	pvda_delete_backend(module);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_univ_tnx(mid_tnx);
	delete_glwe_params(params_glwe);
}
