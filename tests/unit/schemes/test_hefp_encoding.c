#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <sys/types.h>

#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "hefp.h"
#include "hefp_private.h"
#include "maths_structures.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

PvdaParamTest(hefp_internal_1, back_and_forth, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);
	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	PolyUnivRnX* interm_vec  = new_univ_rnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	rnx_random_vec(module, initial_vec, params_glwe);

	hefp_encode_slow_internal(interm_vec, nn, (complex double*)initial_vec);
	hefp_decode_slow_internal((complex double*)final_vec, nn, interm_vec);

	pvda_assert_polynomial_distance(params_glwe, initial_vec, final_vec, 0.0001, 0.0001);

	delete_univ_rnx(initial_vec);
	delete_univ_rnx(final_vec);
	delete_univ_rnx(interm_vec);
	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(hefp_fft, forth_and_back_slow, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);
	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	PolyUnivRnX* interm_vec  = new_univ_rnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	rnx_random_vec(module, initial_vec, params_glwe);

	hefp_encode_internal(module, interm_vec, nn / 2, (complex double*)initial_vec, 0);
	hefp_decode_slow_internal((complex double*)final_vec, nn, interm_vec);

	pvda_assert_polynomial_distance(params_glwe, initial_vec, final_vec, 0.0001, 0.0001);

	delete_univ_rnx(initial_vec);
	delete_univ_rnx(final_vec);
	delete_univ_rnx(interm_vec);
	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(hefp_fft, forth_slow_and_back, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);
	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	PolyUnivRnX* interm_vec  = new_univ_rnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	rnx_random_vec(module, initial_vec, params_glwe);

	hefp_encode_slow_internal(interm_vec, nn, (complex double*)initial_vec);
	hefp_decode_internal(module, (complex double*)final_vec, nn / 2, interm_vec);

	pvda_assert_polynomial_distance(params_glwe, initial_vec, final_vec, 0.0001, 0.0001);

	delete_univ_rnx(initial_vec);
	delete_univ_rnx(final_vec);
	delete_univ_rnx(interm_vec);
	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(hefp_fft, forth_and_back, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);
	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	PolyUnivRnX* interm_vec  = new_univ_rnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	rnx_random_vec(module, initial_vec, params_glwe);

	hefp_encode_internal(module, interm_vec, nn / 2, (complex double*)initial_vec, 0);
	hefp_decode_internal(module, (complex double*)final_vec, nn / 2, interm_vec);

	pvda_assert_polynomial_distance(params_glwe, initial_vec, final_vec, 0.0001, 0.0001);

	delete_univ_rnx(initial_vec);
	delete_univ_rnx(final_vec);
	delete_univ_rnx(interm_vec);
	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(hefp_encoding, forth_and_back, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);
	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	PolyBiv* interm_vec      = new_biv(params_glwe);
	PolyUnivTnX* tmp_tnx     = new_univ_tnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	// (ab)use the fact that tnx and Z mod 2^64 are isomorphic and
	// the memory representation is the same for isomprphic values
	// We need to limit the bits to less than 64 due to the shift we will do
	uniform_pow2_random_pol_znx(module, (PolyUniv*)tmp_tnx, params_glwe->nn, 52);

	// Test for different scaling values
	for (int i = -9; i <= 0; ++i)
	{
		univ_tnx_to_rnx(params_glwe, initial_vec, tmp_tnx);

		hefp_encode(module, params_glwe, interm_vec, nn / 2, i, (complex double*)initial_vec);
		hefp_decode(module, params_glwe, (complex double*)final_vec, nn / 2, i, interm_vec);

		// Compare the two vectors
		pvda_assert_polynomial_distance(params_glwe, initial_vec, final_vec, 0.0001, 0.0001);
	}

	delete_univ_rnx(initial_vec);
	delete_univ_rnx(final_vec);
	delete_biv(interm_vec);
	delete_univ_tnx(tmp_tnx);

	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(hefp_encoding, forth_and_back_not_all_slots, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);
	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	PolyBiv* interm_vec      = new_biv(params_glwe);
	PolyUnivTnX* tmp_tnx     = new_univ_tnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	// (ab)use the fact that tnx and Z mod 2^64 are isomorphic and
	// the memory representation is the same for isomprphic values
	// We need to limit the bits to less than 64 due to the shift we will do
	uniform_pow2_random_pol_znx(module, (PolyUniv*)tmp_tnx, params_glwe->nn, 52);

	// Test for different number of slots, from 2 to nn/2
	for (int i = 2; i < nn; i *= 2)
	{
		univ_tnx_to_rnx(params_glwe, initial_vec, tmp_tnx);

		hefp_encode(module, params_glwe, interm_vec, i, 0, (complex double*)initial_vec);
		hefp_decode(module, params_glwe, (complex double*)final_vec, i, 0, interm_vec);

		// Create a new params object so that we only compare the first "i*2" coefficients,
		// which is actually the first i complex slots
		GLWEParams compar_params = *params_glwe;
		compar_params.nn         = i * 2;
		pvda_assert_polynomial_distance(&compar_params, initial_vec, final_vec, 0.0001, 0.0001);
	}

	delete_univ_rnx(initial_vec);
	delete_univ_rnx(final_vec);
	delete_biv(interm_vec);
	delete_univ_tnx(tmp_tnx);

	DELETE_PVDA_PARAMS_GGSW;
}
