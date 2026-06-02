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

#define MATRIX_COLS  4
#define LOG2_COLS    2
#define MATRIX_ROWS  4

#define NBASE        (1 << 11)
#define KBASE        1
#define KAPPABASE    19
#define LBASE        5
#define NLIMBSBASE   (LBASE * 2)
#define LGADBASE     10
#define LGGBASE      10
#define NGGLIMBSBASE (LGGBASE * 2)

#define SHFT_AMT     51
int onionpir_prepare_query_rows() {}

void print_coefs_glwe(const MODULE* module, const GLWECiphertext* glwe, const GLWESecretKeyPrepared* sk_prep, int n,
                      int shft)
{
	PolyBiv* result_biv     = new_biv_poly(glwe->params);
	PolyUnivTnX* result_tnx = new_univ_tnx(glwe->params);
	glwe_secret_decrypt(module, result_biv, sk_prep, glwe);
	biv_to_univ_tnx(glwe->params, result_tnx, result_biv);

	for (int i = 0; i < n; ++i)
	{
		printf("%lx (%ld)  ", result_tnx[i], ((result_tnx[i] >> (shft - 1)) + 1) >> 1);
	}

	delete_biv(result_biv);
	delete_univ_tnx(result_tnx);
}

int onionpir_fill_bivariate_with_matrix_position(const GLWEParams* params_glwe, PolyBiv* biv, int64_t row,
                                                 int64_t column)
{
	PolyUniv* test = new_univ(params_glwe);
	int64_t rn     = (row + 1l) << SHFT_AMT;
	int64_t cn     = (column + 1l) << SHFT_AMT;
	//cn             = (9l) << SHFT_AMT;
	//cn             = row << SHFT_AMT;

	printf("Filling R %ld C %ld rn %lx cn %lx sr %ld rc %ld    ", row, column, rn, cn, rn >> SHFT_AMT, cn >> SHFT_AMT);

	for (int i = 0; i < NBASE; ++i)
	{
		if (i % 2)
			test[i] = cn;
		else
			test[i] = rn;

		if (i < 4) printf("%ld ", test[i] >> SHFT_AMT);
	}
	printf("\n");
	univ_tnx_to_biv(params_glwe, biv, test, 0);

	delete_univ(test);
	return 0;
}

int main()
{
	double sigma            = ldexp(1.0, 2 - (LBASE)*KAPPABASE);
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* params_glwegad = new_glwegadget_params(params_glwe, KAPPABASE, LGADBASE);
	GGSWParams* params_ggsw          = new_ggsw_params(params_glwe, KBASE, KAPPABASE, NGGLIMBSBASE);

	//Client phase 1
	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);

	uniform_glwe_secret_key(module, sk, 2);
	glwe_sk_prepare(module, sk_prep, sk);

	// Automorphism-expand keys
	GLWEAutomorphismKSKCollection* ksks = new_automorphism_ksk_collection(2 * params_glwe->nn);
	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKSK* ksk = new_automorphism_ksk(params_glwegad);
		prepare_automorphism_key(module, ksk, sk_prep, p);
		glwegadget_ksk_collection_put_key(ksks, ksk, p);
	}

	// GGSW(-s) for gadget to GGSW conversion
	GGSWCiphertextPrep** ggsw_ksks = calloc(KBASE, sizeof(GGSWCiphertextPrep*));
	PolyUniv* neg_sk_i             = new_univ(params_glwe);
	GGSWCiphertext* tmp_ggsw       = new_ggsw(params_ggsw);
	for (uint64_t i = 0; i < KBASE; ++i)
	{
		ggsw_ksks[i] = new_ggsw_prep(params_ggsw);
		for (int p = 0; p < NBASE; ++p)
		{
			neg_sk_i[p] = -glwe_prepared_sk_extract_poly_coefs(sk_prep, i)[p];
		}
		ggsw_secret_encrypt(module, tmp_ggsw, sk_prep, neg_sk_i);
		ggsw_prepare(module, ggsw_ksks[i], tmp_ggsw);
	}
	delete_ggsw(tmp_ggsw);
	free(neg_sk_i);

	GLWECiphertext* row_query = new_glwe(params_glwe);
	GLWECiphertext* col_query = new_glwe(params_glwe);

	PolyUniv* sel_row = new_univ(params_glwe);
	PolyUniv* sel_col = new_univ(params_glwe);
	memset(sel_row, 0, poly_univ_bytes(params_glwe));
	memset(sel_col, 0, poly_univ_bytes(params_glwe));
	sel_col[0] = 1;
	sel_col[1] = 1;
	sel_row[3] = 1;

	glwegadget_packed_secret_encrypt(module, row_query, params_glwegad, sk_prep, sel_row, MATRIX_ROWS);
	glwegadget_packed_secret_encrypt(module, col_query, params_glwegad, sk_prep, sel_col, LOG2_COLS);

	//Server
	struct timespec server_start;
	clock_gettime(CLOCK_REALTIME, &server_start);

	GLWECiphertext* glwe_tree[LOG2_COLS + 1][MATRIX_COLS] = {0};
	uint64_t used_cols                                    = MATRIX_COLS;
	for (int c = 0; c < used_cols; ++c)
	{
		glwe_tree[0][c] = new_glwe(params_glwe);
	}

	GLWEGadgetCiphertextPrep** glwegad_trace = calloc(MATRIX_ROWS, sizeof(GLWEGadgetCiphertextPrep*));
	glwegad_trace[0]                         = new_glwegadget_prep(params_glwegad);
	packed_glwegadget_trace_expand_prepared(module, glwegad_trace, MATRIX_ROWS, row_query, ksks);

	//Half products
	PolyBiv* pos_biv         = new_biv_poly(params_glwe);
	GLWECiphertext* tmp_glwe = new_glwe(params_glwe);
	for (int64_t c = 0; c < MATRIX_COLS; ++c)
	{
		for (int64_t r = 0; r < MATRIX_ROWS; ++r)
		{
			onionpir_fill_bivariate_with_matrix_position(params_glwe, pos_biv, r, c);
			glwegadget_half_prod(module, tmp_glwe, glwegad_trace[r], pos_biv);
			add_glwe(module, glwe_tree[0][c], glwe_tree[0][c], tmp_glwe);
		}
		printf("\n");
	}
	for (int64_t c = 0; c < MATRIX_COLS; ++c)
	{
		printf("Col %d: ", c);
		print_coefs_glwe(module, glwe_tree[0][c], sk_prep, 4, SHFT_AMT);
		printf("\n");
	}

	for (int r = 0; r < MATRIX_ROWS; ++r)
	{
		delete_glwegadget_prep(glwegad_trace[r]);
	}

	GGSWCiphertextPrep* ggsw_trace[LOG2_COLS] = {0};
	ggsw_trace[0]                             = new_ggsw_prep(params_ggsw);

	packed_glwegadget_trace_expand_ggsw_prepared(module, ggsw_trace, LOG2_COLS, ggsw_params_l_tilde_a(params_ggsw),
	                                             col_query, ksks, (const GGSWCiphertextPrep**)ggsw_ksks);

	//CMux tree
	printf("\n");
	for (int l = 0; l < LOG2_COLS; ++l)
	{
		//Init next level
		for (int c = 0; c < (used_cols + 1) / 2; ++c)
		{
			glwe_tree[l + 1][c] = new_glwe(params_glwe);
		}

		//CMux
		int c;
		printf("Lvl %d:\n", l);
		for (c = 0; c < used_cols; c += 2)
		{
			if (c + 1 < used_cols)
				tfhe_cmux(module, glwe_tree[l + 1][c / 2], glwe_tree[l][c], glwe_tree[l][c + 1], ggsw_trace[l], 1);
			else
				tfhe_cmux(module, glwe_tree[l + 1][c / 2], glwe_tree[l][c], glwe_tree[l][c], ggsw_trace[l], 1);
			printf("Col %d: ", c);
			print_coefs_glwe(module, glwe_tree[l + 1][c / 2], sk_prep, 4, SHFT_AMT);
			printf("\n");
		}
		printf("\n");

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
	struct timespec server_end;
	clock_gettime(CLOCK_REALTIME, &server_end);

	printf("Server elapsed time: %ld ms\n",
	       (server_end.tv_sec - server_start.tv_sec) * 1000 + (server_end.tv_nsec - server_start.tv_nsec) / 1000000);

	print_coefs_glwe(module, glwe_tree[LOG2_COLS][0], sk_prep, 4, SHFT_AMT);

	delete_glwe(glwe_tree[LOG2_COLS][0]);

	return 0;
}
