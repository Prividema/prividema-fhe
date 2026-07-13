#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <sys/types.h>

#include "ckks.h"
#include "maths_structures.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

PvdaParamTest(ckks_internal_1, back_and_forth, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);
	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	PolyUnivRnX* interm_vec  = new_univ_rnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	rnx_random_vec(module, initial_vec, params_glwe);

	encode_slow_internal(interm_vec, nn, (complex double*)initial_vec);
	decode_slow_internal((complex double*)final_vec, nn, interm_vec);

	pvda_assert_polynomial_distance(params_glwe, initial_vec, final_vec, 0.0001, 0.0001);

	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(ckks_fft, forth_and_back_slow, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);
	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	PolyUnivRnX* interm_vec  = new_univ_rnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	rnx_random_vec(module, initial_vec, params_glwe);

	encode_internal(module, interm_vec, nn / 2, (complex double*)initial_vec, 0);
	decode_slow_internal((complex double*)final_vec, nn, interm_vec);

	pvda_assert_polynomial_distance(params_glwe, initial_vec, final_vec, 0.0001, 0.0001);

	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(ckks_fft, forth_slow_and_back, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);
	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	PolyUnivRnX* interm_vec  = new_univ_rnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	rnx_random_vec(module, initial_vec, params_glwe);

	encode_slow_internal(interm_vec, nn, (complex double*)initial_vec);
	decode_internal(module, (complex double*)final_vec, nn / 2, interm_vec);

	pvda_assert_polynomial_distance(params_glwe, initial_vec, final_vec, 0.0001, 0.0001);

	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(ckks_fft, forth_and_back, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);
	PolyUnivRnX* initial_vec = new_univ_rnx(params_glwe);
	PolyUnivRnX* final_vec   = new_univ_rnx(params_glwe);
	PolyUnivRnX* interm_vec  = new_univ_rnx(params_glwe);

	uint64_t nn = params_glwe->nn;

	rnx_random_vec(module, initial_vec, params_glwe);

	encode_internal(module, interm_vec, nn / 2, (complex double*)initial_vec, 0);
	decode_internal(module, (complex double*)final_vec, nn / 2, interm_vec);

	pvda_assert_polynomial_distance(params_glwe, initial_vec, final_vec, 0.0001, 0.0001);

	DELETE_PVDA_PARAMS_GGSW;
}
