#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdlib.h>

#include "backend.h"
#include "backend_arithmetic.h"
#include "core/glwe/glwe_params.h"
#include "test_utils.h"

// DFT a 3-polynomial vector, take a view at each offset, iDFT the view and check it matches the input from that offset.
PvdaParamTest(biv_dft_at, view_matches_input_limbs, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	const uint64_t nn   = params_glwe->nn;
	const uint64_t size = 3;

	int64_t* in  = malloc(size * nn * sizeof(int64_t));
	int64_t* out = malloc(size * nn * sizeof(int64_t));
	BE_PolyBivDFT* vec_dft = pvda_new_biv_dft(module, size);
	cr_assert(in != NULL && out != NULL && vec_dft != NULL);

	// Small values so the FFT64 round trip is exact
	for (uint64_t i = 0; i < size * nn; ++i) in[i] = (rand() % 2001) - 1000;

	PolyBiv in_vec = {nn, size, (int64_t)nn, in};
	cr_assert(pvda_biv_to_dft(module, vec_dft, size, &in_vec) == 0);

	for (uint64_t offset = 0; offset < size; ++offset)
	{
		uint64_t view_size = size - offset;

		BE_PolyBivDFT* view             = pvda_biv_dft_at(module, vec_dft, offset);
		const BE_PolyBivDFT* view_const = pvda_biv_dft_at_const(module, vec_dft, offset);
		cr_assert(view == view_const, "const and non-const views differ at offset %lu", offset);

		PolyBiv out_vec = {nn, view_size, (int64_t)nn, out};
		cr_assert(pvda_biv_idft(module, &out_vec, view_const, view_size) == 0);

		for (uint64_t i = 0; i < view_size * nn; ++i)
			cr_assert(eq(i64, out[i], in[offset * nn + i]), "offset %lu, coefficient %lu", offset, i);
	}

	pvda_delete_biv_dft(module, vec_dft);
	free(in);
	free(out);

	DELETE_PVDA_PARAMS_GLWE;
}
