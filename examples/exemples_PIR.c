#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bivariate_polynomial.h"
#include "ggsw_arithmetic.h"
#include "ggsw_ciphertext.h"
#include "ggsw_key.h"
#include "ggsw_params.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "glwegadget_arithmetic.h"
#include "glwegadget_ciphertext.h"
#include "glwegadget_key.h"
#include "maths_structures.h"
#include "schemes/tfhe.h"
#include "spqlios_alias.h"
#include "univariate_polynomial.h"
#include "utils.h"

// #define MATRIX_COLS 16384
// #define LOG2_COLS   14
#define MATRIX_COLS        256
#define LOG2_COLS          8
#define MATRIX_ROWS        256

#define IN_MEMORY_DFT_COLS 64

#define NBASE              (1 << 12)
#define KBASE              1
#define KAPPABASE          18

#define SHFT_AMT           10

#define L_TILDE_Q1         4

#define PRINTPARTIAL       (0)

//Global vairable for debugging
GLWESecretKeyPrepared* dbg_key = NULL;

// Bivariate polynomial print funciton for both debugging and the example
void print_coefs_biv(const PolyBiv* biv, int max_n, int max_l)
{
	for (int l = 0; l < max_l && l < biv->l; ++l)
	{
		printf("Biv lvl %03d: ", l);

		for (int p = 0; p < max_n && p < biv->nn; ++p)
		{
			printf("%03ld ", biv->ptr[l * biv->stride + p]);
		}
		printf("\n");
	}
}

// Function that generates the data for a certain matrix position
// Replace with desired placeholder for the example, or with acual data
// when going to production
int onionpir_fill_bivariate_with_matrix_position(const GLWEParams* params_glwe, PolyBiv* biv, int64_t row,
                                                 int64_t column)
{
	PolyUnivTnX* test = new_univ_tnx(params_glwe);
	uint64_t rn       = (row + 1l) << SHFT_AMT;
	uint64_t cn       = (column + 1l) << SHFT_AMT;

	for (int i = 0; i < NBASE; ++i)
	{
		if (i % 2)
			test[i] = cn;
		else
			test[i] = rn;
	}
	univ_tnx_to_biv(params_glwe, biv, (uint64_t*)test, 0);

	delete_univ_tnx(test);
	return 0;
}

// Creates a matrix-column bivariate polynomial by generating each of its elements and concatenating them
int onionpir_fill_column_with_matrix_position(const GLWEParams* params_glwe, PolyBiv* biv, int64_t column,
                                              int64_t biv_depth, int64_t rows)
{
	assert(biv->l == biv_depth * rows);

	for (int i = 0; i < rows; ++i)
	{
		PolyBiv rowbiv = {biv->nn, biv_depth, biv->nn, biv->ptr + (biv_depth * biv->nn) * i};
		onionpir_fill_bivariate_with_matrix_position(params_glwe, &rowbiv, i, column);
	}
}

// Prepared column variable
// In a real use case, this would not exist and we would instead use a
// memory map of some kind
PolyBivDFT* columns_dft[IN_MEMORY_DFT_COLS + 1] = {0};

// Retrieve a prepared column. This would be simple pointer arithmetic
// to an mmaped file in disk in a real usecase
PolyBivDFT* onionpir_get_prepared_column(int64_t column)
{
	if (column >= IN_MEMORY_DFT_COLS) return columns_dft[0];
	return columns_dft[column + 1];
}

// Function to pre-process some of the rows of the database
// For the example, only the first IN_MEMORY_DFT_COLS are pre-processed for memory limitation reasons
int prepare_column(const MODULE* module, int64_t column, const GLWEParams* db_params,
                   const GLWEGadgetParams* query1_params)
{
	assert(query1_params->l_tilde == glwe_params_l_a(db_params));
	uint64_t total_depth = query1_params->l_tilde * MATRIX_ROWS;
	if (column >= IN_MEMORY_DFT_COLS)
	{
		if (columns_dft[0] == NULL)
		{
			columns_dft[0] = new_biv_dft_custom_params(db_params->nn, total_depth);
			//Already memset to 0 at init time
		}
		return 0;
	}
	PolyBivDFT* pos_biv_dft = new_biv_dft_custom_params(db_params->nn, total_depth);
	PolyBiv* pos_biv        = new_biv_custom_params(db_params->nn, total_depth);
	columns_dft[column + 1] = pos_biv_dft;

	onionpir_fill_column_with_matrix_position(db_params, pos_biv, column, query1_params->l_tilde, MATRIX_ROWS);

	GLWEParams total_params = {db_params->nn, db_params->k, db_params->kappa,
	                           MATRIX_ROWS * db_params->ciphertext_nb_limbs};

	biv_coefs_to_prep(module, &total_params, pos_biv_dft, pos_biv);
	delete_biv(pos_biv);
}

// Performs the server tasks in OnionPIR: receive the packed GLWEGadgets,
// unpack them (expand them), do the Half-product per each column, and finally
// select a column with the CMux tree
int onionpir_server(const MODULE* module, const GGSWParams* ggsw_ksk_params, const GLWEGadgetParams* query1_params,
                    const GLWEParams* db_params, const GLWEParams* aggregation_params,
                    const GLWEAutomorphismKSKCollection* ksks, const GGSWCiphertextPrep** ggsw_ksks,
                    GLWECiphertext* res, const GLWECiphertext* row_query, const GLWECiphertext* col_query)
{
	GLWECiphertext* glwe_tree[LOG2_COLS + 1][MATRIX_COLS] = {0};
	uint64_t used_cols                                    = MATRIX_COLS;
	for (int c = 0; c < used_cols; ++c)
	{
		glwe_tree[0][c] = new_glwe(aggregation_params);
	}

	GLWEGadgetParams* mega_params = new_glwegadget_params(query1_params->params_glwe, query1_params->kappa_tilde,
	                                                      query1_params->l_tilde * MATRIX_ROWS);
	GLWEGadgetCiphertextPrep* glwegad_trace = new_glwegadget_prep(mega_params);

	int st = packed_glwegadget_trace_expand_prepared_single(module, glwegad_trace, query1_params, MATRIX_ROWS,
	                                                        L_TILDE_Q1, row_query, ksks);
	//Half products
	GLWECiphertextDFT* tmp_glwe_dft = new_glwe_dft(aggregation_params);
	struct timespec server_start;
	clock_gettime(CLOCK_REALTIME, &server_start);
	for (int64_t c = 0; c < MATRIX_COLS; ++c)
	{
		glwegadget_half_prod_prepared_to_dft(module, tmp_glwe_dft, glwegad_trace, onionpir_get_prepared_column(c));
		glwe_dft_to_coef(module, glwe_tree[0][c], tmp_glwe_dft);
	}
	struct timespec server_end;
	clock_gettime(CLOCK_REALTIME, &server_end);

	double ms_elapsed =
	    (server_end.tv_sec - server_start.tv_sec) * 1000 + (server_end.tv_nsec - server_start.tv_nsec) / 1000000;
	printf("HP elapsed time: %.2f ms\n", ms_elapsed);

	delete_glwegadget_prep(glwegad_trace);

	GGSWCiphertextPrep* ggsw_trace[LOG2_COLS] = {0};

	packed_glwegadget_trace_expand_ggsw_prepared(module, ggsw_trace, ggsw_ksk_params, LOG2_COLS,
	                                             ggsw_params_l_tilde_a(ggsw_ksk_params), col_query, ksks,
	                                             (const GGSWCiphertextPrep**)ggsw_ksks);

	tfhe_cmux_tree(module, res, (const GLWECiphertext**)&glwe_tree[0], MATRIX_COLS,
	               (const GGSWCiphertextPrep**)ggsw_trace, LOG2_COLS, 1);

	for (int i = 0; i < LOG2_COLS; ++i)
	{
		delete_ggsw_prep(ggsw_trace[i]);
	}
	delete_glwegadget_params(mega_params);
	return 0;
}

// Setup phase for the client in the protocol: secret and evaluation key generation
int onionpir_client_phase0(MODULE* module, GLWESecretKeyPrepared** sk_prep_out, int sk_bits, GLWEParams* sk_params,
                           GLWEAutomorphismKSKCollection** ksks_out, const GLWEGadgetParams* auto_ksk_params,
                           GGSWCiphertextPrep*** ggsw_ksks_out, const GGSWParams* auto_ggsw_params)
{
	int status                          = -1;
	GLWESecretKey* sk                   = alloc_glwe_secret_key(sk_params);
	GLWESecretKeyPrepared* sk_prep      = alloc_glwe_secret_key_prepared(sk_params);
	GLWEAutomorphismKSKCollection* ksks = new_automorphism_ksk_collection(2ul * NBASE);
	GGSWCiphertextPrep** ggsw_ksks      = (GGSWCiphertextPrep**)calloc(KBASE, sizeof(GGSWCiphertextPrep*));
	CHECK_ALLOC(sk, "Secret key allocation failed in phase0 of OnionPIR");
	CHECK_ALLOC(sk_prep, "Prepared secret key allocation failed in phase0 of OnionPIR");
	CHECK_ALLOC(ksks, "Allocation failed in phase0 of onionPIR");
	CHECK_ALLOC(ggsw_ksks, "Allocation failed in phase 0 of onionPIR");

	//Secret key
	*sk_prep_out = sk_prep;
	dbg_key      = sk_prep;
	uniform_glwe_secret_key(module, sk, sk_bits);
	glwe_sk_prepare(module, sk_prep, sk);

	// Automorphism-expand keys
	*ksks_out = ksks;
	for (uint64_t i = 1; (1ULL << i) <= NBASE; ++i)
	{
		int64_t p                = (int64_t)NBASE / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKSK* ksk = new_automorphism_ksk(auto_ksk_params);
		prepare_automorphism_key(module, ksk, sk_prep, (int)p);
		glwegadget_ksk_collection_put_key(ksks, ksk, p);
	}

	// GGSW(-s) for gadget to GGSW conversion
	*ggsw_ksks_out = ggsw_ksks;
	generate_glwegad_to_ggsw_ksk(module, ggsw_ksks, auto_ggsw_params, sk_prep);

	status = 0;
cleanup:
	delete_glwe_secret_key(sk);
	return status;
}

// Initial client phase: generate the row and column packed GLWEGadgets according to the
// desired row and column to select
int onionpir_client_phase1(const MODULE* module, GLWECiphertext** row_query, GLWECiphertext** col_query,
                           const GLWESecretKeyPrepared* sk_prep, int row, int column,
                           const GLWEParams* params_row_query, const GLWEParams* params_col_query,
                           const GLWEGadgetParams* row_query_gad_params, const GLWEGadgetParams* col_query_gad_params)
{
	int status        = -1;
	*row_query        = new_glwe(params_row_query);
	*col_query        = new_glwe(params_col_query);
	PolyUniv* sel_row = new_univ(params_row_query);
	PolyUniv* sel_col = new_univ(params_col_query);

	memset(sel_row, 0, poly_univ_bytes(params_row_query));
	memset(sel_col, 0, poly_univ_bytes(params_col_query));

	uint64_t col_num = column;
	uint64_t row_num = row;

	for (int i = 0; i < LOG2_COLS; ++i)
	{
		sel_col[i] = (int64_t)col_num % 2;
		col_num >>= 1;
	}

	sel_row[row_num] = 1;

	glwegadget_packed_secret_encrypt(module, *row_query, row_query_gad_params, sk_prep, sel_row, MATRIX_ROWS);
	glwegadget_packed_secret_encrypt(module, *col_query, col_query_gad_params, sk_prep, sel_col, LOG2_COLS);

	status = 0;
cleanup:
	delete_univ(sel_row);
	delete_univ(sel_col);
	return status;
}

int main()
{
	double sigma8  = ldexp(1.0, 2 - 8 * KAPPABASE);
	double sigma5  = ldexp(1.0, 2 - 5 * KAPPABASE);
	double sigma4  = ldexp(1.0, 2 - 4 * KAPPABASE);
	double sigma6  = ldexp(1.0, 2 - 6 * KAPPABASE);
	MODULE* module = pvda_new_module_info(NBASE);

	// Automorphims keys
	double ksk_sigma = sigma8;
	GLWEParams* params_glwe_autokey =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, 16, ksk_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* auto_ksk_params = new_glwegadget_params(params_glwe_autokey, KAPPABASE, 8);

	//GGSW keys
	double ggsw_ksk_sigma = sigma8;
	GLWEParams* params_ggsw_change_key =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, 16, ggsw_ksk_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GGSWParams* auto_ggsw_params = new_ggsw_params(params_glwe_autokey, KBASE, KAPPABASE, 16);

	//Input queries
	//col query
	double row_sigma             = sigma8;
	GLWEParams* params_row_query = new_glwe_params(NBASE, KBASE, KAPPABASE, 16, row_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* row_query_gad_params = new_glwegadget_params(params_row_query, KAPPABASE, L_TILDE_Q1);
	//col query
	double col_sigma             = sigma8;
	GLWEParams* params_col_query = new_glwe_params(NBASE, KBASE, KAPPABASE, 16, col_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* col_query_gad_params = new_glwegadget_params(params_col_query, KAPPABASE, 8);

	//Expanded queries
	// row expansion
	double row_exp_sigma = sigma6;
	GLWEParams* row_exp_params =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, 12, row_exp_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* row_exp_gad_params = new_glwegadget_params(row_exp_params, KAPPABASE, L_TILDE_Q1);

	// col expansion
	double col_exp_sigma = sigma6;
	GLWEParams* col_exp_params =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, 12, col_exp_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GGSWParams* col_exp_ggsw_params = new_ggsw_params(col_exp_params, KBASE, KAPPABASE, 10);

	double col_sum_sigma = sigma6;
	GLWEParams* col_sum_params =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, 12, col_sum_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEParams* final_params = new_glwe_params(NBASE, KBASE, KAPPABASE, 10, sigma5, NOISE_UNIFORM_POWER_OF_TWO);

	GLWEParams* db_params = new_glwe_params(NBASE, KBASE, KAPPABASE, 8, 0, NOISE_UNIFORM_POWER_OF_TWO);

	//Client phase 0

	GLWESecretKeyPrepared* sk_prep;
	GLWEAutomorphismKSKCollection* ksks;
	GGSWCiphertextPrep** ggsw_ksks;

	onionpir_client_phase0(module, &sk_prep, 2, final_params, &ksks, auto_ksk_params, &ggsw_ksks, auto_ggsw_params);

	GLWECiphertext* row_query;
	GLWECiphertext* col_query;

	onionpir_client_phase1(module, &row_query, &col_query, sk_prep, 3, 10, params_row_query, params_col_query,
	                       row_query_gad_params, col_query_gad_params);

	GLWECiphertext* res = new_glwe(final_params);

	for (int c = 0; c <= IN_MEMORY_DFT_COLS; ++c)
	{
		prepare_column(module, c, db_params, row_exp_gad_params);
	}

	struct timespec server_start;
	clock_gettime(CLOCK_REALTIME, &server_start);

	onionpir_server(module, auto_ggsw_params, row_exp_gad_params, db_params, col_sum_params, ksks,
	                (const GGSWCiphertextPrep**)ggsw_ksks, res, row_query, col_query);

	struct timespec server_end;
	clock_gettime(CLOCK_REALTIME, &server_end);

	double ms_elapsed =
	    (server_end.tv_sec - server_start.tv_sec) * 1000 + (server_end.tv_nsec - server_start.tv_nsec) / 1000000;
	printf("Server elapsed time: %.2f ms\n", ms_elapsed);

	print_coefs_glwe(module, res, sk_prep, 4, SHFT_AMT);

	double throughput_bits_sec = (64.0 - SHFT_AMT) * NBASE * MATRIX_ROWS * MATRIX_COLS * 1000 / ms_elapsed;
	printf("Server throughput: %.2f MiB/s\n", throughput_bits_sec / 8 / 1024 / 1024);
	fflush(stdout);

	delete_glwe(res);

	return 0;
}
