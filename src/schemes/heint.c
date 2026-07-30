#include "heint.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>

#include "backend.h"
#include "backend_private.h"
#include "bivariate_polynomial.h"
#include "math_utils.h"
#include "montgomery_arith.h"
#include "pvda_ffts.h"
#include "univariate_polynomial.h"
#include "utils.h"

int internal_slow_intt_heint(const PvdaBackend* backend, uint64_t* out_int, uint64_t* in, uint64_t t)
{
	int status  = -1;
	uint64_t nn = pvda_module_extract_nn(backend);

	uint64_t t_tild        = montgomery_tild_32bit(t);
	uint64_t r2modt        = montgomery_r2modq_32bit(t);
	uint64_t log2n         = next_pow2_log(nn + 1);
	uint64_t twice_nn_mask = (1ull << log2n) - 1;
	uint64_t n_inv         = mod_inv(nn, t);

	uint64_t n_inv_m = montgomery_encode_32bit(n_inv, t, t_tild, r2modt);

	NTTRoot* root_table    = get_ntt_table(backend->pvda_fft_data, t);
	uint64_t* root_table_m = malloc((2 * nn + 1) * sizeof(uint64_t));
	uint64_t* in_m         = malloc(nn * sizeof(uint64_t));
	CHECK_ALLOC(root_table, "root table not found");
	CHECK_ALLOC(in_m, "failed malloc in heint encoding");
	CHECK_ALLOC(root_table_m, "failed malloc in heint encoding");

	// Prepare root table
	for (size_t i = 0; i < 2 * nn + 1; ++i) root_table_m[i] = montgomery_encode_32bit(root_table[i], t, t_tild, r2modt);

	for (size_t i = 0; i < nn; ++i)
	{
		in_m[i] = montgomery_encode_32bit(in[i], t, t_tild, r2modt);
	}

	for (size_t i = 0; i < nn; ++i)
	{
		// 32-bit montgomery arithmetic
		uint32_t sum = 0;

		for (size_t j = 0; j < nn / 2; ++j)
		{
			uint64_t root_num = backend->pvda_fft_data->rotation_group[j] * i;
			uint32_t p1_red   = montgomery_mult_32bit(in_m[j], root_table_m[(-root_num) & twice_nn_mask], t, t_tild);

			uint32_t p2_red =
			    montgomery_mult_32bit(in_m[j + nn / 2], root_table_m[root_num & twice_nn_mask], t, t_tild);

			sum += p1_red;
			sum += p2_red;
		}

		uint64_t sum_n_m = (uint64_t)sum * n_inv_m;
		uint64_t sum_red = montgomery_red_32bit(sum_n_m, t, t_tild);
		uint64_t sum_dec = montgomery_decode_32bit(sum_red, t, t_tild);
		out_int[i]       = sum_dec;
	}
	status = 0;
cleanup:
	free(in_m);
	free(root_table_m);
	return status;
}

int internal_slow_ntt_heint(const PvdaBackend* backend, uint64_t* out_int, uint64_t* in, uint64_t t)
{
	int status = -1;

	uint64_t nn            = pvda_module_extract_nn(backend);
	uint64_t t_tild        = montgomery_tild_32bit(t);
	uint64_t r2modt        = montgomery_r2modq_32bit(t);
	uint64_t log2n         = next_pow2_log(nn + 1);
	uint64_t twice_nn_mask = (1ull << log2n) - 1;

	NTTRoot* root_table    = get_ntt_table(backend->pvda_fft_data, t);
	uint64_t* root_table_m = malloc((2 * nn + 1) * sizeof(uint64_t));
	uint64_t* in_m         = malloc(nn * sizeof(uint64_t));
	CHECK_ALLOC(root_table, "root table not found");
	CHECK_ALLOC(in_m, "failed malloc in heint encoding");
	CHECK_ALLOC(root_table_m, "failed malloc in heint encoding");

	// Prepare root table
	for (size_t i = 0; i < 2 * nn + 1; ++i) root_table_m[i] = montgomery_encode_32bit(root_table[i], t, t_tild, r2modt);

	for (size_t i = 0; i < nn; ++i)
	{
		in_m[i] = montgomery_encode_32bit(in[i], t, t_tild, r2modt);
	}

	for (size_t i = 0; i < nn / 2; ++i)
	{
		// 32-bit montgomery arithmetic
		uint32_t sum_1 = 0;
		uint32_t sum_2 = 0;

		for (size_t j = 0; j < nn; ++j)
		{
			uint64_t root_num = backend->pvda_fft_data->rotation_group[i] * j;
			uint64_t air1     = montgomery_mult_32bit(in_m[j], root_table_m[(root_num)&twice_nn_mask], t, t_tild);
			uint64_t air2     = montgomery_mult_32bit(in_m[j], root_table_m[(-root_num) & twice_nn_mask], t, t_tild);

			// Radix R (32-bit) addition
			sum_1 += air1;
			sum_2 += air2;
		}

		out_int[i]          = montgomery_decode_32bit(sum_1, t, t_tild);
		out_int[i + nn / 2] = montgomery_decode_32bit(sum_2, t, t_tild);
	}
	status = 0;
cleanup:
	free(in_m);
	free(root_table_m);
	return status;
}

int heint_encode(const PvdaBackend* backend, const GLWEParams* params, PolyBiv* out, uint64_t slots, int64_t t,
                 uint64_t* in)
{
	int status  = -1;
	uint64_t nn = pvda_module_extract_nn(backend);

	assert(t <= (1ll << 31));

	NTTRoot* root_table    = get_ntt_table(backend->pvda_fft_data, t);
	uint64_t t_tild        = montgomery_tild_32bit(t);
	uint64_t r2modt        = montgomery_r2modq_32bit(t);
	uint64_t* root_table_m = malloc((2 * nn + 1) * sizeof(uint64_t));
	uint64_t* out_int      = malloc(nn * sizeof(uint64_t));
	uint64_t* out_tnx      = malloc(nn * sizeof(uint64_t));
	double* out_rnx        = malloc(nn * sizeof(double));

	CHECK_ALLOC(out_int, "malloc failed in heint encoding");
	CHECK_ALLOC(root_table_m, "failed malloc in heint encoding");
	CHECK_ALLOC(root_table, "root table not found");
	CHECK_ALLOC(out_tnx, "malloc failed in heint encoding");
	CHECK_ALLOC(out_rnx, "failed malloc in heint encoding");

	if (slots != nn)
	{
		RAISE_ERROR("Unsupported (for now) number of slots != N");
	}

	// Prepare roots in montgomery form
	for (size_t i = 0; i < 2 * nn + 1; ++i)
	{
		root_table_m[i] = montgomery_encode_32bit(root_table[i], t, t_tild, r2modt);
	}

	// Perform iNTT
	internal_slow_intt_heint(backend, out_int, in, t);

	// Scale output (BFV style)
	double scale_fact = 1.0 / t;
	for (size_t i = 0; i < nn; ++i) out_rnx[i] = scale_fact * out_int[i];

	// Since all positions can share the same scale factor,
	// go rnx -> tnx -> biv since it is faster than the
	// direct rnx -> biv conversion (which supports different scales)
	univ_rnx_to_tnx(params, out_tnx, out_rnx);
	univ_tnx_to_biv(params, out, out_tnx, 0);

	status = 0;
cleanup:
	free(root_table_m);
	free(out_int);
	free(out_tnx);
	free(out_rnx);

	return status;
}

int heint_decode(const PvdaBackend* backend, const GLWEParams* params, uint64_t* out, uint64_t slots, int64_t t,
                 PolyBiv* in)
{
	int status  = -1;
	uint64_t nn = pvda_module_extract_nn(backend);

	assert(t <= (1ll << 31));

	uint64_t* tmp_tnx = malloc(nn * sizeof(uint64_t));
	double* in_rnx    = malloc(nn * sizeof(double));

	CHECK_ALLOC(tmp_tnx, "malloc failed in heint encoding");
	CHECK_ALLOC(in_rnx, "failed malloc in heint encoding");

	if (slots != nn) RAISE_ERROR("Unsupported (for now) number of slots != N");

	biv_to_univ_rnx(params, in_rnx, in, 0);

	// Use the [0,1] torus representation
	for (size_t i = 0; i < nn; ++i) in_rnx[i] = in_rnx[i] < 0 ? 1 + in_rnx[i] : in_rnx[i];

	// Scale and round
	for (size_t i = 0; i < nn; ++i) tmp_tnx[i] = llround(t * in_rnx[i]);

	internal_slow_ntt_heint(backend, out, tmp_tnx, t);
	status = 0;
cleanup:
	free(tmp_tnx);
	free(in_rnx);

	return status;
}
