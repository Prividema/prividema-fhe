#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <stdint.h>
#include <sys/types.h>

#include "backend.h"
#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "heint.h"
#include "maths_structures.h"
#include "pvda_ffts.h"
#include "rng.h"
#include "univariate_polynomial.h"

Test(heint_helpers, mont_pow)
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
}
