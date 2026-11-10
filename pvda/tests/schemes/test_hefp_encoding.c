#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <sys/types.h>

#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "hefp.h"
#include "hefp_private.h"
#include "random_polynomial.h"
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
	// the memory representation is the same for isomorphic values
	// We need to limit the bits to less than 64 due to the shift we will do
	univ_znx_random_uniform_pow2(module, (int64_t*)tmp_tnx, params_glwe->nn, 52);

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
	// the memory representation is the same for isomorphic values
	// We need to limit the bits to less than 64 due to the shift we will do
	univ_znx_random_uniform_pow2(module, (int64_t*)tmp_tnx, params_glwe->nn, 52);

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

// =============================================
// |                                           |
// |              hefp_decode_fft              |
// |                                           |
// =============================================

#include <complex.h>
#include <math.h>

/** @brief Runs hefp_decode_fft on the coefficients in (packed as re = in[i], im = in[i + N/2]) for ring size nn. */
static void decode_fft_vs_slow(uint64_t nn, const double* in)
{
	PvdaBackend* module   = pvda_new_spqlios_backend(nn);
	uint64_t slots        = nn / 2;
	double complex* fast  = malloc(slots * sizeof(double complex));
	double complex* slow  = malloc(slots * sizeof(double complex));
	for (uint64_t i = 0; i < slots; ++i) fast[i] = CMPLX(in[i], in[i + slots]);

	cr_assert(eq(int, hefp_decode_fft(module, fast, slots), 0));
	hefp_decode_slow_internal(slow, nn, (double*)in);
	for (uint64_t j = 0; j < slots; ++j)
		cr_assert(cabs(fast[j] - slow[j]) < 1e-9 * (double)nn, "nn %lu slot %lu: %g%+gi vs %g%+gi", nn, j,
		          creal(fast[j]), cimag(fast[j]), creal(slow[j]), cimag(slow[j]));

	free(slow);
	free(fast);
	pvda_delete_backend(module);
}

Test(hefp_decode_fft, matches_slow_decoding)
{
	static const uint64_t sizes[] = {4, 8, 64, 1024};  // 4: the smallest case, slots = 2
	for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s)
	{
		uint64_t nn = sizes[s];
		double* in  = malloc(nn * sizeof(double));
		for (uint64_t i = 0; i < nn; ++i) in[i] = sin(1.0 + 3.7 * (double)i) * 0.5;
		decode_fft_vs_slow(nn, in);
		free(in);
	}
}

Test(hefp_decode_fft, known_polynomials)
{
	const uint64_t nn   = 16, slots = nn / 2;
	PvdaBackend* module = pvda_new_spqlios_backend(nn);
	double complex v[8];

	// The constant polynomial 1 evaluates to 1 in every slot
	for (uint64_t i = 0; i < slots; ++i) v[i] = i == 0 ? 1 : 0;
	cr_assert(eq(int, hefp_decode_fft(module, v, slots), 0));
	for (uint64_t j = 0; j < slots; ++j) cr_assert(cabs(v[j] - 1) < 1e-12, "slot %lu", j);

	// The zero polynomial gives 0
	for (uint64_t i = 0; i < slots; ++i) v[i] = 0;
	cr_assert(eq(int, hefp_decode_fft(module, v, slots), 0));
	for (uint64_t j = 0; j < slots; ++j) cr_assert(eq(dbl, cabs(v[j]), 0.0));

	// X evaluates to the primitive 2N-th roots exp(i pi 5^j / N), one per slot
	for (uint64_t i = 0; i < slots; ++i) v[i] = i == 1 ? 1 : 0;
	cr_assert(eq(int, hefp_decode_fft(module, v, slots), 0));
	uint64_t pw = 1;
	for (uint64_t j = 0; j < slots; ++j)
	{
		double complex root = cexp(I * M_PI * (double)pw / (double)nn);
		cr_assert(cabs(v[j] - root) < 1e-12, "slot %lu", j);
		pw = (pw * 5) % (2 * nn);
	}

	pvda_delete_backend(module);
}
