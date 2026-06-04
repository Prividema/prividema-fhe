#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bivariate_polynomial.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe_arithmetic.h"
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
#define MATRIX_COLS  1024
#define LOG2_COLS    10
#define MATRIX_ROWS  1024

#define NBASE        (1 << 12)
#define KBASE        1
#define KAPPABASE    16

#define SHFT_AMT     16

#define L_TILDE_Q1   4

#define PRINTPARTIAL (0)

GLWESecretKeyPrepared* dbg_key = NULL;

int onionpir_prepare_query_rows() {}

void print_coefs_glwe(const MODULE* module, const GLWECiphertext* glwe, const GLWESecretKeyPrepared* sk_prep, int n,
                      int shft)
{
	PolyBiv* result_biv     = new_biv_poly(glwe->params);
	PolyUnivTnX* result_tnx = new_univ_tnx(glwe->params);
	//normalize_glwe(module, glwe, glwe);
	glwe_secret_decrypt(module, result_biv, sk_prep, glwe);
	biv_to_univ_tnx(glwe->params, result_tnx, result_biv);

	for (int i = 0; i < n; ++i)
	{
		if (shft)
			printf("%lx (%ld)  ", result_tnx[i], ((result_tnx[i] >> (shft - 1)) + 1) >> 1);
		else
			printf("%ld ", result_tnx[i]);
	}

	delete_biv(result_biv);
	delete_univ_tnx(result_tnx);
}

void print_coefs_gad(const MODULE* module, const GLWEGadgetCiphertext* glwe_gad, const GLWESecretKeyPrepared* sk_prep,
                     int n)
{
	for (int l = 0; l < glwe_gad->params->l_tilde; ++l)
	{
		GLWECiphertext limb = {glwe_gad->params->params_glwe, glwegadget_extract_bivglwe(glwe_gad, l + 1)};
		printf("gad l %03d: ", l);
		print_coefs_glwe(module, &limb, sk_prep, n, 0);
		printf("\n");
	}
}

void print_coefs_biv(const PolyBiv* biv, int max_n, int max_l)
{
	for (int l = 0; l < max_l && l < biv->l; ++l)
	{
		printf("biv lvl %03d: ", l);

		for (int p = 0; p < max_n && p < biv->nn; ++p)
		{
			printf("%03ld", biv->ptr[l * biv->stride + p]);
		}
		printf("\n");
	}
}

int onionpir_fill_bivariate_with_matrix_position(const GLWEParams* params_glwe, PolyBiv* biv, int64_t row,
                                                 int64_t column)
{
	PolyUniv* test = new_univ(params_glwe);
	int64_t rn     = (row + 1l) << SHFT_AMT;
	int64_t cn     = (column + 1l) << SHFT_AMT;

	for (int i = 0; i < NBASE; ++i)
	{
		if (i % 2)
			test[i] = cn;
		else
			test[i] = rn;
	}
	univ_tnx_to_biv(params_glwe, biv, test, 0);

	delete_univ(test);
	return 0;
}

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
	//packed_glwegadget_trace_expand_prepared(module, glwegad_trace, MATRIX_ROWS, 3, row_query, ksks);

	if (PRINTPARTIAL)
	{
		GLWEGadgetCiphertext** glwegad_trace_unprep = calloc(MATRIX_ROWS, sizeof(GLWEGadgetCiphertext*));
		for (int i = 0; i < MATRIX_ROWS; ++i) glwegad_trace_unprep[i] = new_glwegadget(query1_params);
		st = packed_glwegadget_trace_expand(module, glwegad_trace_unprep, MATRIX_ROWS, 4, row_query, ksks);

		printf("Unprep : ");
		print_coefs_gad(module, glwegad_trace_unprep[18], dbg_key, 4 * MATRIX_ROWS);
		printf("\n");

		for (int i = 0; i < MATRIX_ROWS; ++i) delete_glwegadget(glwegad_trace_unprep[i]);
	}

	// printf("st=%d\n", st);
	//Half products
	GLWECiphertext* tmp_glwe = new_glwe(aggregation_params);

	PolyBiv* pos_biv = new_biv_poly_custom_l(db_params, query1_params->l_tilde * MATRIX_ROWS);
	for (int64_t c = 0; c < MATRIX_COLS; ++c)
	{
		//
		onionpir_fill_column_with_matrix_position(db_params, pos_biv, c, query1_params->l_tilde, MATRIX_ROWS);
		glwegadget_half_prod(module, glwe_tree[0][c], glwegad_trace, pos_biv);

		if (PRINTPARTIAL)
		{
			printf("Column %d: ", c);
			print_coefs_glwe(module, glwe_tree[0][c], dbg_key, 4, SHFT_AMT);
			printf("\n");
		}
	}
	delete_glwegadget_prep(glwegad_trace);

	GGSWCiphertextPrep* ggsw_trace[LOG2_COLS] = {0};
	ggsw_trace[0]                             = new_ggsw_prep(ggsw_ksk_params);

	packed_glwegadget_trace_expand_ggsw_prepared(module, ggsw_trace, LOG2_COLS, ggsw_params_l_tilde_a(ggsw_ksk_params),
	                                             col_query, ksks, (const GGSWCiphertextPrep**)ggsw_ksks);

	//CMux tree
	if (PRINTPARTIAL) printf("\n");
	for (int l = 0; l < LOG2_COLS; ++l)
	{
		//Init next level
		if (l == LOG2_COLS - 1)
			glwe_tree[l + 1][0] = res;
		else
			for (int c = 0; c < (used_cols + 1) / 2; ++c)
			{
				glwe_tree[l + 1][c] = new_glwe(aggregation_params);
			}

		//CMux
		int c;
		if (PRINTPARTIAL) printf("Lvl %d:\n", l);
		for (c = 0; c < used_cols; c += 2)
		{
			if (c + 1 < used_cols)
				tfhe_cmux(module, glwe_tree[l + 1][c / 2], glwe_tree[l][c], glwe_tree[l][c + 1], ggsw_trace[l], 1);
			else
			{
				glwe_copy(glwe_tree[l + 1][c / 2], glwe_tree[l][c]);
			}
			if (PRINTPARTIAL)
			{
				printf("Column %d: ", c / 2);
				print_coefs_glwe(module, glwe_tree[l + 1][c / 2], dbg_key, 4, SHFT_AMT);
				printf("\n");
			}
		}
		if (PRINTPARTIAL) printf("\n");

		for (int c = 0; c < used_cols; ++c)
		{
			delete_glwe(glwe_tree[l][c]);
		}
		used_cols = (used_cols + 1) / 2;
	}

	for (int i = 0; i < LOG2_COLS; ++i)
	{
		delete_ggsw_prep(ggsw_trace[i]);
	}
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

	//GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_UNIFORM_POWER_OF_TWO);
	// GLWEGadgetParams* params_glwegad = new_glwegadget_params(params_glwe, KAPPABASE, LGADBASE);
	// GGSWParams* params_ggsw          = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NGGLIMBSBASE);

	GLWEParams* db_params = new_glwe_params(NBASE, KBASE, KAPPABASE, 6, 0, NOISE_UNIFORM_POWER_OF_TWO);
	//Client phase 1
	GLWESecretKey* sk              = alloc_glwe_secret_key(final_params);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(final_params);
	dbg_key                        = sk_prep;

	uniform_glwe_secret_key(module, sk, 2);
	glwe_sk_prepare(module, sk_prep, sk);

	// Automorphism-expand keys
	GLWEAutomorphismKSKCollection* ksks = new_automorphism_ksk_collection(2 * NBASE);
	for (uint64_t i = 1; (1ULL << i) <= NBASE; ++i)
	{
		int64_t p                = (int64_t)NBASE / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKSK* ksk = new_automorphism_ksk(auto_ksk_params);
		prepare_automorphism_key(module, ksk, sk_prep, p);
		glwegadget_ksk_collection_put_key(ksks, ksk, p);
	}

	// GGSW(-s) for gadget to GGSW conversion
	GGSWCiphertextPrep** ggsw_ksks = calloc(KBASE, sizeof(GGSWCiphertextPrep*));
	PolyUniv* neg_sk_i             = new_univ(params_ggsw_change_key);
	GGSWCiphertext* tmp_ggsw       = new_ggsw(auto_ggsw_params);
	for (uint64_t i = 0; i < KBASE; ++i)
	{
		ggsw_ksks[i] = new_ggsw_prep(auto_ggsw_params);
		for (int p = 0; p < NBASE; ++p)
		{
			neg_sk_i[p] = -glwe_prepared_sk_extract_poly_coefs(sk_prep, i)[p];
		}
		ggsw_secret_encrypt(module, tmp_ggsw, sk_prep, neg_sk_i);
		ggsw_prepare(module, ggsw_ksks[i], tmp_ggsw);
	}
	delete_ggsw(tmp_ggsw);
	free(neg_sk_i);

	GLWECiphertext* row_query = new_glwe(params_row_query);
	GLWECiphertext* col_query = new_glwe(params_col_query);

	PolyUniv* sel_row = new_univ(params_row_query);
	PolyUniv* sel_col = new_univ(params_col_query);
	memset(sel_row, 0, poly_univ_bytes(params_row_query));
	memset(sel_col, 0, poly_univ_bytes(params_col_query));
	sel_col[0]  = 1;
	sel_col[1]  = 1;
	sel_col[5]  = 1;
	sel_row[17] = 1;

	glwegadget_packed_secret_encrypt(module, row_query, row_query_gad_params, sk_prep, sel_row, MATRIX_ROWS);
	glwegadget_packed_secret_encrypt(module, col_query, col_query_gad_params, sk_prep, sel_col, LOG2_COLS);

	printf("Row query: ");
	print_coefs_glwe(module, row_query, sk_prep, MATRIX_ROWS * 4, 0);
	printf("\n");

	GLWECiphertext* res = new_glwe(final_params);
	//Server

	struct timespec server_start;
	clock_gettime(CLOCK_REALTIME, &server_start);

	onionpir_server(module, auto_ggsw_params, row_exp_gad_params, db_params, col_sum_params,

	                ksks, (const GGSWCiphertextPrep**)ggsw_ksks, res, row_query, col_query);

	struct timespec server_end;
	clock_gettime(CLOCK_REALTIME, &server_end);

	double ms_elapsed =
	    (server_end.tv_sec - server_start.tv_sec) * 1000 + (server_end.tv_nsec - server_start.tv_nsec) / 1000000;
	printf("Server elapsed time: %.2f ms\n", ms_elapsed);

	print_coefs_glwe(module, res, sk_prep, 4, SHFT_AMT);

	double throughput_bits_sec = (64.0 - SHFT_AMT) * NBASE * MATRIX_ROWS * MATRIX_COLS * 1000 / ms_elapsed;
	printf("Server throughput: %.2f MB/s\n", throughput_bits_sec / 8 / 1024 / 2024);
	fflush(stdout);

	delete_glwe(res);

	return 0;
}
