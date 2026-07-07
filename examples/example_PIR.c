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

/*****************************************************************************
 * OionPIR example using Half-products:
 *
 * This file contains an example implementation of OnionPIR
 * (or possibly more closely OnionPIRv2) that leverages the double decomposition
 * technique and the base-2K representation together with the half-product.
 *
 * See:
 * OnionPIR: https://eprint.iacr.org/2021/1081
 * OnionPIRv2: https://eprint.iacr.org/2025/1142
 *
 *****************************************************************************/

// Matrix dimenstions (number of columns, rows and logarithm of the columns)
#define MATRIX_COLS 1024
#define LOG2_COLS   10
#define MATRIX_ROWS 1024

// How many columns to actually keep in memory
// In a real application, all of them would be either in disk or memory
// However, for development and proof-of-concept, we only use the first
// IN_MEMORY_DFT_COLS and for the rest we reuse a 0 column
// Otherwise, we would quickly run out of memory

#define IN_MEMORY_DFT_COLS 32

// GLWE(and Gadget) parameters
#define NBASE      4096
#define KBASE      1
#define KAPPABASE  18
#define L_TILDE_Q1 4

// Elements in the database are technically in the torus.
// For simplicity, in this example the data we would like to put is
// some number in ZnX mod 2^64 or equivalently TnX
// However, due to the noisy nature of FHE, the result we get at the end
// will have some noise. To avoid that, we would want the data in the DB to
// be BFV-like, that is, have some margin in the LSBs where the noise will
// go and which we then can just discard and round to get an exact value
//
// SHFT_AMT is this margin, in bits. More succintly, in this OnionPIR proof-of-concept
// we store DB_BITS bits of useful data per coefficient, and everything after that
// is margin for the noise
//
// Even more succintly, the data in the database belongs in ZnX mod 2^(DB_BITS)
#define DB_BITS  54ll
#define SHFT_AMT (64 - DB_BITS)

// Size of the square submatrix for which we print the DB contents at the beggining
// of the exampole
#define DB_PRINT_SUBMAT 4
#define DB_PRINT_COEFS  4

// Function that generates the data for a certain matrix position
// Replace with desired placeholder for the example, or with acual data
// when going to production
// At present, it places a polynomial with alternaticng column number and row number coefficents
int onionpir_fill_bivariate_with_matrix_position(const GLWEParams* params_glwe, PolyBiv* biv, int64_t row,
                                                 int64_t column)
{
	PolyUnivTnX* test = new_univ_tnx(params_glwe);
	uint64_t rn       = (row + 1l) << SHFT_AMT;
	uint64_t cn       = (column + 1l) << SHFT_AMT;

	for (int i = 0; i < NBASE; ++i)
	{
		//Sample pattern of alternating column and row number in the polynomial coefficents.
		if (i % 2)
			test[i] = cn;
		else
			test[i] = rn;
	}
	if (row < DB_PRINT_SUBMAT && column < DB_PRINT_SUBMAT)
	{
		printf("DB pos %ld %ld contents:", row + 1, column + 1);
		for (int i = 0; i < DB_PRINT_COEFS; ++i) printf("   %ld (%ld)", test[i], test[i] >> SHFT_AMT);
		printf("...\n");
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

// This stores the database columns after being pre-processed ("prepared") for faster operation
// In a real use case, this would not exist as a variable in memory.
// Instead, we would use a memory mapping and it would reside in disk
PolyBivDFT* columns_dft[IN_MEMORY_DFT_COLS + 1] = {0};

// Retrieve a prepared column. This would be simple pointer arithmetic
// to an mmaped file in disk in a real use case
PolyBivDFT* onionpir_get_prepared_column(int64_t column)
{
	if (column >= IN_MEMORY_DFT_COLS) return columns_dft[0];
	return columns_dft[column + 1];
}

// Function to pre-process some of the columns of the database
//
// The database is stored in a "prepared" DFT format that improves the performance of queries
//
// For the example, only the first IN_MEMORY_DFT_COLS are pre-processed for memory limitation reasons
// Column 0 is stored in columns_dft[1], and so on until columns_dft[IN_MEMORY_DFT_COLS]
// columns_dft[0] is used to store a 0-valued column placeholder used for columns that are
// greater than IN_MEMORY_DFT_COLS
int prepare_column(const MODULE* module, int64_t column, const GLWEParams* db_params,
                   const GLWEGadgetParams* query1_params)
{
	int status = -1;
	assert(query1_params->l_tilde == glwe_params_l_a(db_params));
	uint64_t total_depth = query1_params->l_tilde * MATRIX_ROWS;
	if (column >= IN_MEMORY_DFT_COLS)
	{
		// If greater than the prepared columns, it should use a 0-valued column
		// If it has not yet been set, do it, otherwise just return
		if (columns_dft[0] == NULL)
		{
			columns_dft[0] = new_biv_dft_custom_params(db_params->nn, total_depth);
			// We do not need to do memset to 0 at this point since new_biv_dft_custom_params already
			// sets the data to 0
			// Unsure if this will change in the future for performance reasons
			// In any case, writing to a DFT polynomial directly might not be supported in platforms
			// where DFT elements, for example, reside in a GPU, etc.
			if (columns_dft[0] == 0)
			{
				log_message(LOG_ERROR, "prepare column 0-column allocation failed");
				return -1;
			}
		}
		return 0;
	}

	// Custom-depth bivariate polynomials, since each "bivariate polynomial" actually corresponds
	// to the concatenation of all the bivariate polynomials in the column
	PolyBivPrep* pos_biv_dft = new_biv_dft_custom_params(db_params->nn, total_depth);
	PolyBiv* pos_biv         = new_biv_custom_params(db_params->nn, total_depth);
	CHECK_ALLOC(pos_biv_dft, "column 'bivariate' polynomial (in DFT) allocation failed");
	CHECK_ALLOC(pos_biv, "column 'bivariate' polynomial allocation failed");

	columns_dft[column + 1] = pos_biv_dft;

	// Fill a whole column with the data generated or retrieved by the generator functions
	CHECK_CALL(
	    onionpir_fill_column_with_matrix_position(db_params, pos_biv, column, query1_params->l_tilde, MATRIX_ROWS),
	    "Filling the onionPIR column failed");

	// This GLWEParams object is a dummy parameter set that is used  only for biv_coefs_to_prep to be able to work on a
	// "bivariate_polynomial" of depth total_dept of depth total_depth.
	//
	GLWEParams total_params = {db_params->nn, db_params->k, db_params->kappa,
	                           MATRIX_ROWS * db_params->ciphertext_nb_limbs};

	//Prepare the column
	CHECK_CALL(biv_coefs_to_prep(module, &total_params, pos_biv_dft, pos_biv),
	           "Column preparation (preprocessing) failed");
	status = 0;
cleanup:
	delete_biv(pos_biv);
	return status;
}

// Performs the server tasks in OnionPIR: receive the packed GLWEGadgets,
// unpack them (expand them), do the Half-product per each column, and finally
// select a column with the CMux tree
int onionpir_server(const MODULE* module, const GGSWParams* ggsw_ksk_params, const GLWEGadgetParams* query1_params,
                    const GLWEParams* db_params, const GLWEParams* aggregation_params,
                    const GLWEAutomorphismKeyCollection* auto_key_collection, const GGSWCiphertextPrep** ggsw_ksks,
                    GLWECiphertext* res, const GLWECiphertext* row_query, const GLWECiphertext* col_query)
{
	int status = -1;

	// New parameters used for the concatenation of GLWEGadgets x Concatenation of rows in column optimisation
	// We have to have a GLWEGadget of depth MATRIX_ROWS times the original depth
	GLWEGadgetParams* mega_params = new_glwegadget_params(query1_params->params_glwe, query1_params->kappa_tilde,
	                                                      query1_params->l_tilde * MATRIX_ROWS);

	// The glwegad_trace is the set of GLWEGadgets, one per row, that encrypt either a 0 or 1 to select a row
	// (note that theoretically, they could encode other values if we wanted a linear combination too)
	// The optimisation used consists in representing said set as a single giant GLWEGadget, since
	// using that we can convert MATRIX_ROWS vector x matrix products plus (MATRIX_ROWS - 1) additions in the DFT domain
	// into a single vector x matrix product, with the dimension of both the matrix and vector being MATRIX_ROWS times as large
	// The theoretical complexity is the same, but reducing it to one operation instead of many can allow for easier optimisation.
	//
	// For example, due to spqlios not having (at present) implementations for addition in the DFT space, we would have had to
	// convert the results of each half-product (vector-matrix product) into coefficient space for aggregation,
	// requiring us to do MATRIX_ROWS DFT-to-coefficient domain conversions instead of a single one at the end
	GLWEGadgetCiphertextPrep* glwegad_trace = new_glwegadget_prep(mega_params);

	// Temporary DFT glwe for the result for each column
	GLWECiphertextDFT* tmp_glwe_dft = new_glwe_dft(aggregation_params);

	CHECK_ALLOC_LABEL(mega_params, "Parameter allocation failed in onionpir server", cleanup1);
	CHECK_ALLOC_LABEL(glwegad_trace, "GLWEGadget trace allocation failed in onionpir server", cleanup1);
	CHECK_ALLOC_LABEL(tmp_glwe_dft, "Temporary GLWE DFT failed in onionpir server", cleanup1);

	// Expand the row selection packed glwegadget into a single, large depth, GLWEGadget
	int st = packed_glwegadget_trace_expand_prepared_single(module, glwegad_trace, query1_params, MATRIX_ROWS,
	                                                        L_TILDE_Q1, row_query, auto_key_collection);
	CHECK_CALL_LABEL(st, "Row trace expansion failed in onionpir server", cleanup1);

	// Half products
	// Due to the optimisation, we only have one per column

	GLWECiphertext* glwe_tree_first_level[MATRIX_COLS] = {0};
	for (int c = 0; c < MATRIX_COLS; ++c)
	{
		glwe_tree_first_level[c] = new_glwe(aggregation_params);
		CHECK_ALLOC_LABEL(glwe_tree_first_level[c], "GLWE allocation failed in onionpir server", cleanup2);
	}
	// We record and report the time spent in the column half-products, since for non-small matrices it is the main performance bottleneck
	struct timespec server_start;
	clock_gettime(CLOCK_REALTIME, &server_start);

	for (int64_t c = 0; c < MATRIX_COLS; ++c)
	{
		// We do just one half-product per column, and then convert it to coefficient-space for later selection with the CMux tree
		st = glwegadget_half_prod_prepared_to_dft(module, tmp_glwe_dft, glwegad_trace, onionpir_get_prepared_column(c));
		CHECK_CALL_LABEL(st, "Trace by column half-product failed in onionpir server", cleanup2);
		st = glwe_dft_to_coef(module, glwe_tree_first_level[c], tmp_glwe_dft);
		CHECK_CALL_LABEL(st, "DFT to coefficient form failed after column half-product in onionpir server", cleanup2);
	}

	// Delete the row selection gadget to lower peak memory consumption
	delete_glwegadget_prep(glwegad_trace);
	glwegad_trace = NULL;

	struct timespec server_end;
	clock_gettime(CLOCK_REALTIME, &server_end);

	double ms_elapsed =
	    (server_end.tv_sec - server_start.tv_sec) * 1000 + (server_end.tv_nsec - server_start.tv_nsec) / 1000000;
	printf("HP elapsed time: %.2f ms\n", ms_elapsed);

	// Column selection by a Mux tree. Thus, we have LOG2_COLS selection signals (one signal per level)
	GGSWCiphertextPrep* ggsw_trace[LOG2_COLS] = {0};

	// Expand the column selection packed GLWEGadget into one GGSW per level (LOG2_COLS as many of them)
	// Note that we do not initialize the ggsw_trace contents.
	// The function has support for doing it automatically, which allows it to reduce the memory footprint
	// (see comment in packed_glwegadget_trace_expand_ggsw_prepared)
	//
	// Notable observations:
	// - auto_key_collection are the keys used for automorphisms
	// - ggsw_ksks are encryptions of the coefficients of -sk, used to (internally) convert the GLWEGadgets that result of the query
	//   expansion into GGSWs
	st = packed_glwegadget_trace_expand_ggsw_prepared(module, ggsw_trace, ggsw_ksk_params, LOG2_COLS,
	                                                  ggsw_params_l_tilde_a(ggsw_ksk_params), col_query,
	                                                  auto_key_collection, (const GGSWCiphertextPrep**)ggsw_ksks);
	CHECK_CALL_LABEL(st, "Column selection trace expansion failed in onionpir server", cleanup3);

	// CMux selection tree (that is, arbitrary-size CMux)
	// This function can deallocate its input, which we do for memory efficiency reasons
	st = tfhe_cmux_tree(module, res, (const GLWECiphertext**)&glwe_tree_first_level, MATRIX_COLS,
	                    (const GGSWCiphertextPrep**)ggsw_trace, LOG2_COLS, 1);
	CHECK_CALL_LABEL(st, "CMux tree failed in onionpir server", cleanup3);

	status = 0;
cleanup3:
	for (int i = 0; i < LOG2_COLS; ++i)
	{
		delete_ggsw_prep(ggsw_trace[i]);
	}
cleanup2:
	for (int c = 0; c < MATRIX_COLS; ++c)
	{
		delete_glwe(glwe_tree_first_level[c]);
	}
cleanup1:
	delete_glwegadget_prep(glwegad_trace);
	delete_glwe_dft(tmp_glwe_dft);
	delete_glwegadget_params(mega_params);
	return status;
}

// Setup phase for the client in the protocol: secret and evaluation key generation
int onionpir_client_phase0(MODULE* module, GLWESecretKeyPrepared** sk_prep_out, int sk_bits, GLWEParams* sk_params,
                           GLWEAutomorphismKeyCollection** auto_key_collection_out,
                           const GLWEGadgetParams* auto_ksk_params, GGSWCiphertextPrep*** ggsw_ksks_out,
                           const GGSWParams* auto_ggsw_params)
{
	int status                                         = -1;
	GLWESecretKey* sk                                  = alloc_glwe_secret_key(sk_params);
	GLWESecretKeyPrepared* sk_prep                     = alloc_glwe_secret_key_prepared(sk_params);
	GLWEAutomorphismKeyCollection* auto_key_collection = new_automorphism_key_collection(2ul * NBASE);
	GGSWCiphertextPrep** ggsw_ksks = (GGSWCiphertextPrep**)calloc(KBASE, sizeof(GGSWCiphertextPrep*));
	CHECK_ALLOC(sk, "Secret key allocation failed in phase 0 of OnionPIR");
	CHECK_ALLOC(sk_prep, "Prepared secret key allocation failed in phase 0 of OnionPIR");
	CHECK_ALLOC(auto_key_collection, "Allocation failed in phase 0 of onionPIR");
	CHECK_ALLOC(ggsw_ksks, "Allocation failed in phase 0 of onionPIR");

	*sk_prep_out             = sk_prep;
	*auto_key_collection_out = auto_key_collection;
	*ggsw_ksks_out           = ggsw_ksks;

	// Secret key
	CHECK_CALL(uniform_glwe_secret_key(module, sk, sk_bits), "GLWE secret key generation failed in onionpir phase 0");
	glwe_sk_prepare(module, sk_prep, sk);

	// Automorphism-expand keys
	for (uint64_t i = 1; (1ULL << i) <= NBASE; ++i)
	{
		int64_t p                     = (int64_t)NBASE / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* auto_key = new_automorphism_key(auto_ksk_params);
		CHECK_ALLOC(auto_key, "Automorphism key allocation failed in onionpir phase 0");

		int st = compute_automorphism_key(module, auto_key, sk_prep, (int)p);
		glwegadget_key_collection_put_key(auto_key_collection, auto_key, p);
		CHECK_CALL(st, "Automorphism key computation failed in onionpir phase 0");
	}

	// GGSW(-s) for gadget to GGSW conversion
	CHECK_CALL(generate_glwegad_to_ggsw_ksk(module, ggsw_ksks, auto_ggsw_params, sk_prep),
	           "GLWEGadget to GGSW KSK generation failed in oniopir phase 0");

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

	CHECK_ALLOC(row_query, "GLWE allocation failed in oniopir phase 1");
	CHECK_ALLOC(col_query, "GLWE allocation failed in oniopir phase 1");
	CHECK_ALLOC(sel_row, "Univariate polynomial allocation failed in oniopir phase 1");
	CHECK_ALLOC(sel_col, "Univariate polynomial allocation failed in oniopir phase 1");

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

	CHECK_CALL(
	    glwegadget_packed_secret_encrypt(module, *row_query, row_query_gad_params, sk_prep, sel_row, MATRIX_ROWS),
	    "Row query packing and encryption failed in onionpir phase 1");
	CHECK_CALL(glwegadget_packed_secret_encrypt(module, *col_query, col_query_gad_params, sk_prep, sel_col, LOG2_COLS),
	           "Column query packing and encryption failed in onionpir phase 1");

	status = 0;
cleanup:
	delete_univ(sel_row);
	delete_univ(sel_col);
	return status;
}

int main(int argc, char* argv[])
{
	if (argc != 1 && argc != 3)
	{
		printf(
		    "Usage: example_PIR\n"
		    "       example_PIR row column\n");
		exit(1);
	}

	printf("General parameters: N = %d, k = %d, K = %d\n", NBASE, KBASE, KAPPABASE);
	printf("\n");

	// We set a minimum of 2 bits of noise standard deviation
	// This sets the stdev to be at least 4 if, compared to
	// the standard 3.2.
	// (when we interpret the torus discretized in 2^-L = 2 ^ -(l*K) to the isomorphic Z_q[x] / (x^N + 1) with q = 2^L
	// to be able to use standard security estimations)
	// Due to using the NOISE_UNIFORM_POWER_OF_TWO
	// generator, this will select the closest power of 2 that has a uniform distribution
	// at least as secure as the standard one with the parameter, since
	// Standard: 3.2*sqrt(3) = 5.54 ==> 8
	// 2 bits:   4 * sqrt(3) = 6.9 ==> 8
	//
	// In other words: these values of sigma all correspond to the same security level
	// We need different values since we change the torus precision (similar to modulo switching)
	double sigma8  = ldexp(1.0, 2 - 8 * KAPPABASE);
	double sigma5  = ldexp(1.0, 2 - 5 * KAPPABASE);
	double sigma6  = ldexp(1.0, 2 - 6 * KAPPABASE);
	MODULE* module = pvda_new_module_info(NBASE);

	// Automorphims keys
	double auto_key_sigma = sigma8;
	GLWEParams* params_glwe_autokey =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, 16, auto_key_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* auto_key_params = new_glwegadget_params(params_glwe_autokey, KAPPABASE, 8);

	//GGSW keys
	double ggsw_ksk_sigma = sigma8;
	GLWEParams* params_ggsw_change_key =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, 16, ggsw_ksk_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GGSWParams* auto_ggsw_params = new_ggsw_params(params_glwe_autokey, KBASE, KAPPABASE, 16);

	//Input queries

	//row query (Input GLWEGadget and packed GLWE form)
	double row_sigma             = sigma8;
	GLWEParams* params_row_query = new_glwe_params(NBASE, KBASE, KAPPABASE, 16, row_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* row_query_gad_params = new_glwegadget_params(params_row_query, KAPPABASE, L_TILDE_Q1);
	//col query (Input GLWEGadget and packed GLWE form)
	double col_sigma             = sigma8;
	GLWEParams* params_col_query = new_glwe_params(NBASE, KBASE, KAPPABASE, 16, col_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEGadgetParams* col_query_gad_params = new_glwegadget_params(params_col_query, KAPPABASE, 8);

	//Expanded queries (parameters for the unpacked form)
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

	// Parameters used for the aggregation and/or final result for each column after the half-external products
	double col_sum_sigma = sigma6;
	GLWEParams* col_sum_params =
	    new_glwe_params(NBASE, KBASE, KAPPABASE, 12, col_sum_sigma, NOISE_UNIFORM_POWER_OF_TWO);
	GLWEParams* final_params = new_glwe_params(NBASE, KBASE, KAPPABASE, 10, sigma5, NOISE_UNIFORM_POWER_OF_TWO);

	GLWEParams* db_params = new_glwe_params(NBASE, KBASE, KAPPABASE, 8, 0, NOISE_UNIFORM_POWER_OF_TWO);

	//Client phase 0: secret and evaluation key generation
	GLWESecretKeyPrepared* sk_prep;
	GLWEAutomorphismKeyCollection* auto_key_collection;
	GGSWCiphertextPrep** ggsw_ksks;

	int st = onionpir_client_phase0(module, &sk_prep, 2, final_params, &auto_key_collection, auto_key_params,
	                                &ggsw_ksks, auto_ggsw_params);

	CHECK_CALL(st, "OnionPIR phase 0 (key generation) failed");

	uint64_t db_size_bits = DB_BITS * NBASE * MATRIX_ROWS * MATRIX_COLS;

	printf("Database geometry: %d rows x %d columns matrix\n", MATRIX_ROWS, MATRIX_COLS);
	printf("Element size (cleartext): %.1f KiB\n", (double)(DB_BITS * NBASE) / 1024);
	printf("Database size (cleartext): %.1f MiB\n", (double)db_size_bits / 8 / 1024 / 1024);

	uint64_t prep_db_size_bytes = MATRIX_COLS * MATRIX_ROWS * poly_biv_bytes(db_params);

	printf("Preprocessed DB size (plaintext in prepared form): %.1f MiB\n", (double)prep_db_size_bytes / 1024 / 1024);
	printf(
	    "Note that for this example, only the first %d columns of the DB are being stored in memory, and all others "
	    "are set to 0\n",
	    IN_MEMORY_DFT_COLS);

	printf("\n");

	printf("Row query size: %.1f KiB\n", (double)glwe_params_bytes(params_row_query) / 1024);
	printf("Column query size: %.1f KiB\n", (double)glwe_params_bytes(params_col_query) / 1024);
	printf("Response size: %.1f KiB\n", (double)glwe_params_bytes(final_params) / 1024);

	printf("\n");

	// Server pre-processing: generate the database columns
	GLWECiphertext* res = new_glwe(final_params);
	printf("Filling database...\n");
	for (int c = 0; c <= IN_MEMORY_DFT_COLS; ++c)
	{
		st = prepare_column(module, c, db_params, row_exp_gad_params);
		CHECK_CALL(st, "Preprocessing of a column of the OnionPIR database failed");
	}
	printf("...\nThe full %d x %d matrix has been filled. Only a portion has been printed for readability\n",
	       MATRIX_ROWS, MATRIX_COLS);

	GLWECiphertext* row_query;
	GLWECiphertext* col_query;
	while (1)
	{
		//Client phase 1: row and column query generation

		int row, col;

		if (argc == 3)
		{
			row = atoi(argv[1]);
			col = atoi(argv[2]);
		}
		else
		{
			printf("Input the row and column to select, space-separated: ");
			int st = scanf("%d %d", &row, &col);
			if (st != 2)
			{
				printf("Wrong selection format\n");
				exit(1);
			}
		}

		printf("Selecting row %d, column %d\n", row, col);
		--row;
		--col;

		if (col >= IN_MEMORY_DFT_COLS)
		{
			printf("Note: selected a column that has blank data, expect the output to be 0\n");
		}

		st = onionpir_client_phase1(module, &row_query, &col_query, sk_prep, row, col, params_row_query,
		                            params_col_query, row_query_gad_params, col_query_gad_params);
		CHECK_CALL(st, "OnionPIR phase 1 (query generation from row and column numbers) failed");

		// Server online phase: receive column and row queries, together with evaluation keys (auto_key_collection and ggsw_ksks)
		// Use OnionPIR to retrieve a single database position in the form of a GLWE

		struct timespec server_start;
		clock_gettime(CLOCK_REALTIME, &server_start);

		st = onionpir_server(module, auto_ggsw_params, row_exp_gad_params, db_params, col_sum_params,
		                     auto_key_collection, (const GGSWCiphertextPrep**)ggsw_ksks, res, row_query, col_query);
		CHECK_CALL(st, "OnionPIR server phase failed");

		struct timespec server_end;
		clock_gettime(CLOCK_REALTIME, &server_end);

		double ms_elapsed =
		    (server_end.tv_sec - server_start.tv_sec) * 1000 + (server_end.tv_nsec - server_start.tv_nsec) / 1000000;

		//Client phase 2: receive and decrypt result (done as part of print_coefs_glwe)

		printf("Result     encoded (decoded):\n");
		print_coefs_glwe(module, res, sk_prep, 4, SHFT_AMT);

		// Compute the online throughput, in terms of bits of useful database data (size of the unprepared database)
		// A througput of 500MB/s would mean that a query to a database that stores 1GB of data, which might be using,
		// for example, 4GB due to the format used in OnionPIR, would take 2s to return a result
		//
		// Creating the row and column query and decrypting the result is NOT part of this throughput
		double throughput_bits_sec = db_size_bits * 1000 / ms_elapsed;
		printf("Server elapsed time: %.2f ms\n", ms_elapsed);
		printf("Server throughput (cleartext DB MiB per second): %.2f MiB/s\n", throughput_bits_sec / 8 / 1024 / 1024);
		fflush(stdout);

		delete_glwe(row_query);
		row_query = NULL;
		delete_glwe(col_query);
		col_query = NULL;

		if (argc == 3) break;
		printf("\n");
	}

	delete_glwe(res);
	delete_automorphism_key_collection(auto_key_collection, 1);
	for (int i = 0; i < KBASE; ++i)
	{
		delete_ggsw_prep(ggsw_ksks[i]);
	}
	free(ggsw_ksks);

	delete_glwe_secret_key_prepared(sk_prep);

	pvda_delete_module_info(module);

	delete_glwe_params(params_glwe_autokey);
	delete_glwegadget_params(auto_key_params);

	delete_glwe_params(params_ggsw_change_key);
	delete_ggsw_params(auto_ggsw_params);

	delete_glwe_params(params_row_query);
	delete_glwegadget_params(row_query_gad_params);
	delete_glwe_params(params_col_query);
	delete_glwegadget_params(col_query_gad_params);

	delete_glwe_params(row_exp_params);
	delete_glwegadget_params(row_exp_gad_params);
	delete_glwe_params(col_exp_params);
	delete_ggsw_params(col_exp_ggsw_params);

	delete_glwe_params(col_sum_params);
	delete_glwe_params(final_params);
	delete_glwe_params(db_params);

	return 0;

	// In this case, cleanup is the error path and we simply exit the program
cleanup:
	exit(-1);
}
