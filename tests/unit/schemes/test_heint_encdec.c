#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <stdint.h>
#include <sys/types.h>

#include "backend_private.h"
#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "heint.h"
#include "maths_structures.h"
#include "pvda_ffts.h"
#include "rng.h"
#include "univariate_polynomial.h"

Test(heint_enc, sanity1)
{
	uint64_t nn             = 8;
	PvdaBackend* backend    = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 4, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 17;  // 2*nn + 1

	PolyUnivTnX* in_tnx = new_univ_tnx(params_glwe);
	PolyBiv* biv_out    = new_biv(params_glwe);

	generate_ntt_table(backend, t);
	uniform_random_pol_znx(backend, in_tnx, nn, 0, t - 1);

	int st = heint_encode(backend, params_glwe, biv_out, nn, t, in_tnx);

	cr_assert(eq(i32, st, 0));
	pvda_delete_backend(backend);
	delete_glwe_params(params_glwe);
	delete_univ_tnx(in_tnx);
	delete_biv(biv_out);
}

Test(heint_enc, backforth)
{
	uint64_t nn             = 8;
	PvdaBackend* backend    = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 6, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 17;  // 2*nn + 1

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);
	PolyBiv* biv_out     = new_biv(params_glwe);

	generate_ntt_table(backend, t);
	uniform_random_pol_znx(backend, in_tnx, nn, 0, t - 1);

	int st  = heint_encode(backend, params_glwe, biv_out, nn, t, in_tnx);
	int st2 = heint_decode(backend, params_glwe, out_tnx, nn, t, biv_out);

	cr_assert(eq(i32, st, 0));
	cr_assert(eq(i32, st2, 0));

	for (int i = 0; i < nn; ++i)
	{
		cr_assert(eq(i64, in_tnx[i], out_tnx[i]));
	}

	pvda_delete_backend(backend);
	delete_glwe_params(params_glwe);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_biv(biv_out);
}
Test(heint_enc, backforth3)
{
	uint64_t nn             = 1024;
	PvdaBackend* backend    = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 6, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 4085761;

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);
	PolyBiv* biv_out     = new_biv(params_glwe);

	generate_ntt_table(backend, t);
	uniform_random_pol_znx(backend, in_tnx, nn, 0, t - 1);

	int st  = heint_encode(backend, params_glwe, biv_out, nn, t, in_tnx);
	int st2 = heint_decode(backend, params_glwe, out_tnx, nn, t, biv_out);

	cr_assert(eq(i32, st, 0));
	cr_assert(eq(i32, st2, 0));

	for (int i = 0; i < nn; ++i)
	{
		cr_assert(eq(i64, in_tnx[i], out_tnx[i]));
	}

	pvda_delete_backend(backend);
	delete_glwe_params(params_glwe);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_biv(biv_out);
}

Test(heint_enc, backforth2)
{
	uint64_t nn             = 1024;
	PvdaBackend* backend    = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 6, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 12289;  // 6*2*nn + 1

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);
	PolyBiv* biv_out     = new_biv(params_glwe);

	generate_ntt_table(backend, t);
	uniform_random_pol_znx(backend, in_tnx, nn, 0, t - 1);

	int st  = heint_encode(backend, params_glwe, biv_out, nn, t, in_tnx);
	int st2 = heint_decode(backend, params_glwe, out_tnx, nn, t, biv_out);

	cr_assert(eq(i32, st, 0));
	cr_assert(eq(i32, st2, 0));

	for (int i = 0; i < nn; ++i)
	{
		cr_assert(eq(i64, in_tnx[i], out_tnx[i]));
	}

	pvda_delete_backend(backend);
	delete_glwe_params(params_glwe);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_biv(biv_out);
}

Test(heint_internals, sanity1)
{
	uint64_t nn             = 8;
	PvdaBackend* backend    = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 4, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 17;  // 2*nn + 1

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);

	generate_ntt_table(backend, t);
	uniform_random_pol_znx(backend, in_tnx, nn, 0, t - 1);

	int st = internal_slow_intt_heint(backend, out_tnx, in_tnx, t);

	cr_assert(eq(i32, st, 0));
	pvda_delete_backend(backend);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_glwe_params(params_glwe);
}

Test(heint_internals, sanity2)
{
	uint64_t nn             = 8;
	PvdaBackend* backend    = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 4, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 17;  // 2*nn + 1

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);

	generate_ntt_table(backend, t);
	uniform_random_pol_znx(backend, in_tnx, nn, 0, t - 1);

	int st = internal_slow_ntt_heint(backend, out_tnx, in_tnx, t);

	cr_assert(eq(i32, st, 0));
	pvda_delete_backend(backend);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_glwe_params(params_glwe);
}

Test(heint_internals, backforth)
{
	uint64_t nn             = 8;
	PvdaBackend* backend    = pvda_new_spqlios_backend(nn);
	GLWEParams* params_glwe = new_glwe_params(nn, 1, 19, 4, 0.0001, NOISE_UNIFORM_POWER_OF_TWO);
	uint64_t t              = 17;  // 2*nn + 1

	PolyUnivTnX* in_tnx  = new_univ_tnx(params_glwe);
	PolyUnivTnX* mid_tnx = new_univ_tnx(params_glwe);
	PolyUnivTnX* out_tnx = new_univ_tnx(params_glwe);

	generate_ntt_table(backend, t);
	uniform_random_pol_znx(backend, in_tnx, nn, 0, t - 1);

	int st  = internal_slow_intt_heint(backend, mid_tnx, in_tnx, t);
	int st2 = internal_slow_ntt_heint(backend, out_tnx, mid_tnx, t);
	cr_assert(eq(i32, st, 0));
	cr_assert(eq(i32, st2, 0));

	for (int i = 0; i < nn; ++i)
	{
		cr_assert(eq(i64, in_tnx[i], out_tnx[i]));
	}

	pvda_delete_backend(backend);
	delete_univ_tnx(in_tnx);
	delete_univ_tnx(out_tnx);
	delete_univ_tnx(mid_tnx);
	delete_glwe_params(params_glwe);
}
