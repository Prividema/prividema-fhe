#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdint.h>

#include "bivariate_polynomial.h"
#include "ggsw_ciphertext.h"
#include "rng.h"
#include "univariate_polynomial.h"

#define NBASE            1024
#define KBASE            1
#define KAPPABASE        4
#define NLIMBSBASE       (KBASE + 1) * 2
#define LBASE            NLIMBSBASE / (KBASE + 1)
#define SIGMABASE        1e-7

#define K_TILDEBASE      1
#define KAPPA_TILDEBASE  4
#define NLIMBS_TILDEBASE (K_TILDEBASE + 1) * 1
#define L_TILDEBASE      NLIMBS_TILDEBASE / (K_TILDEBASE + 1)

//! COMMON PART (begin)

// Test ggsw_size
Test(ggsw_size, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Asserts ggsw_size returns NLIMBS_TILDEBASE * NLIMBSBASE
	cr_assert(eq(i64, ggsw_size(params_ggsw), NLIMBS_TILDEBASE * NLIMBSBASE));

	// Clean up
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
}

// Test ggsw_bytes
Test(ggsw_bytes, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Asserts ggsw_bytes returns NLIMBS_TILDEBASE * NLIMBSBASE * NBASE * sizeof(int64_t)
	cr_assert(eq(i64, ggsw_bytes(params_ggsw), NLIMBS_TILDEBASE * NLIMBSBASE * NBASE * sizeof(int64_t)));

	// Clean up
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
}

//! bivGGSW Part (begin)

// Test ggsw_coef_number
Test(ggsw_coef_number, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Asserts ggsw_coeff_number returns NLIMBS_TILDEBASE * NLIMBSBASE * NBASE
	cr_assert(eq(i64, ggsw_coef_number(params_ggsw), NLIMBS_TILDEBASE * NLIMBSBASE * NBASE));

	// Clean up
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
}

// Test new_ggsw
Test(new_ggsw, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	GGSWCiphertext* ggsw = new_ggsw(params_ggsw);

	// Asserts new_ggsw worked
	cr_assert(eq(int, ggsw->mat != NULL, 1));
	cr_assert(eq(int, ggsw->params != NULL, 1));

	// Clean up
	delete_ggsw(ggsw);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
}

// Test ggsw_Sj_Yti
Test(ggsw_Sj_Yti, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	GGSWCiphertext* ggsw = new_ggsw(params_ggsw);

	// asserts ggsw_Sj_Yti returns the right pointer
	for (uint64_t i = 1; i < ggsw_num_pggsw(params_ggsw); i++)
		for (uint64_t j = 0; j < K_TILDEBASE + 1; j++)
		{
			VecBiv* ct_mat_ij = ggsw_retrieve_bivglwe(ggsw, j, i);

			// Modify the two firsts coefficients of biGLWE(-m * sk_j / 2^{kappa_tilde * i}).
			ct_mat_ij[0] = 1;
			ct_mat_ij[1] = 2;

			cr_assert(eq(i64, ggsw->mat[(i - 1) * (K_TILDEBASE + 1) * NLIMBSBASE * NBASE + j * NLIMBSBASE * NBASE], 1));
			cr_assert(
			    eq(i64, ggsw->mat[(i - 1) * (K_TILDEBASE + 1) * NLIMBSBASE * NBASE + j * NLIMBSBASE * NBASE + 1], 2));
		}

	// Clean up
	delete_ggsw(ggsw);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
}

// Test normalize_ggsw
Test(normalize_ggsw, basic)
{
	// Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	GGSWCiphertext* ggsw = new_ggsw(params_ggsw);
	GGSWCiphertext* res  = new_ggsw(params_ggsw);

	// Normalize the ggsw ciphertext
	int status = normalize_ggsw(module, res, ggsw);

	// Asserts normalize_ggsw succeed
	cr_assert(eq(int, status, 0), "normalize_ggsw failed");

	// Clean up
	delete_ggsw(ggsw);
	delete_ggsw(res);
	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
}

Test(add_ggsw, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	GGSWCiphertext* ggsw_lhs     = new_ggsw(params_ggsw);
	GGSWCiphertext* ggsw_rhs     = new_ggsw(params_ggsw);
	GGSWCiphertext* sum_computed = new_ggsw(params_ggsw);

	// Draws uniformly the bivGGSW ciphertexts
	uniform_random_vec(NBASE, ggsw_lhs->mat, ggsw_size(params_ggsw), NBASE, KAPPABASE - 1);
	uniform_random_vec(NBASE, ggsw_rhs->mat, ggsw_size(params_ggsw), NBASE, KAPPABASE - 1);

	// Computes ggsw_lhs + ggsw_rhs
	add_ggsw(sum_computed, ggsw_lhs, ggsw_rhs);

	uint64_t nb_rows = sum_computed->params->n_limbs_tilde;
	uint64_t nb_cols = sum_computed->params->params_glwe->n_limbs;
	uint64_t N       = sum_computed->params->params_glwe->nn;

	// Asserts sum_computed = ggsw_lhs + ggsw_rhs
	for (uint64_t idx = 0; idx < N * nb_cols * nb_rows; ++idx)
		cr_assert(eq(int, sum_computed->mat[idx], ggsw_lhs->mat[idx] + ggsw_rhs->mat[idx]),
		          "add_biv_ggsw mismatch at index %" PRId64 ": %" PRId64 " + %" PRId64 " = %" PRId64 ", got %" PRId64,
		          (long long)idx, ggsw_lhs->mat[idx], ggsw_rhs->mat[idx], ggsw_lhs->mat[idx] + ggsw_rhs->mat[idx],
		          sum_computed->mat[idx]);

	// Clean up
	delete_ggsw(ggsw_lhs);
	delete_ggsw(ggsw_rhs);
	delete_ggsw(sum_computed);
	delete_ggsw_params(params_ggsw);
	delete_glwe_params(params_glwe);
}

Test(const_mult_ggsw, without_normalization)
{
	// Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	PolyUniv* u                      = new_univ(params_glwe);
	PolyUnivDFT* u_dft               = new_univ_dft(module);
	GGSWCiphertext* ggsw             = new_ggsw(params_ggsw);
	GGSWCiphertext* product_computed = new_ggsw(params_ggsw);
	PolyUniv* prod_expected          = new_univ(params_glwe);

	// Draws uniformly the Zn[X] constant
	uniform_random_vec(NBASE, u, 1, NBASE, KAPPABASE - 1);

	// Draws uniformly the bivGGSW ciphertext
	uniform_random_vec(NBASE, ggsw->mat, ggsw_size(params_ggsw), NBASE, KAPPABASE - 1);

	// Computes u in the DFT domain
	univ_coefs_to_dft(module, u_dft, u);

	// Computes u * ggsw
	const_mult_ggsw(module, product_computed, ggsw, u_dft);

	// Asserts product_computed = u * ggsw
	for (uint64_t ii = 1; ii <= ggsw_num_pggsw(params_ggsw); ii++)
		for (uint64_t jj = 0; jj < K_TILDEBASE + 1; jj++)
		{
			VecBiv* ct_mat_ii_jj  = ggsw_retrieve_bivglwe(ggsw, jj, ii);
			VecBiv* res_mat_ii_jj = ggsw_retrieve_bivglwe(product_computed, jj, ii);
			for (uint64_t j = 0; j < KBASE + 1; j++)
				for (uint64_t i = LBASE; i >= 1; i--)
				{
					PolyUniv* ct_ij = ct_mat_ii_jj + (i - 1) * (KBASE + 1) * NBASE + j * NBASE;
					pvda_znx_product(module, prod_expected, u, ct_ij);
					for (uint64_t p = 0; p < NBASE; p++)
					{
						cr_assert(
						    eq(i64, res_mat_ii_jj[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p], prod_expected[p]));
					}
				}
		}

	// Clean up
	delete_univ(u);
	delete_univ(prod_expected);
	delete_univ_dft(u_dft);
	delete_ggsw(ggsw);
	delete_ggsw(product_computed);
	delete_ggsw_params(params_ggsw);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

//! bivGGSW IN DFT SPACE Part (begin)

// Test ggsw_coef_number
Test(ggsw_coef_number_dft, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Asserts ggsw_coef_number_dft returns NLIMBS_TILDEBASE * NLIMBSBASE * NBASE / 2
	cr_assert(eq(i64, ggsw_coef_number_dft(params_ggsw), NLIMBS_TILDEBASE * NLIMBSBASE * NBASE / 2));

	// Clean up
	delete_ggsw_params(params_ggsw);
	delete_glwe_params(params_glwe);
}

// Test new_ggsw
Test(new_ggsw_dft, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	GGSWCiphertextDFT* ggsw_dft = new_ggsw_dft(params_ggsw);

	// Asserts new_ggsw_dft worked
	cr_assert(eq(int, ggsw_dft->mat != NULL, 1));
	cr_assert(eq(int, ggsw_dft->params != NULL, 1));

	// Clean up
	delete_ggsw_dft(ggsw_dft);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
}

// Test ggsw_Sj_Yti_dft
Test(ggsw_Sj_Yti_dft, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	GGSWCiphertextDFT* ggsw_dft = new_ggsw_dft(params_ggsw);

	// Asserts ggsw_Sj_Yti_dft returns the right pointer
	for (uint64_t i = 1; i < ggsw_num_pggsw(params_ggsw); i++)
		for (uint64_t j = 0; j < K_TILDEBASE + 1; j++)
		{
			VecBivDFT* ct_mat_ij = ggsw_retrieve_bivglwe_dft(ggsw_dft, j, i);

			// Modify the two firsts coefficients of biGLWE(-m * sk_j / 2^{kappa_tilde * i}).
			ct_mat_ij[0] = 0.1;
			ct_mat_ij[1] = 0.2;

			cr_assert(
			    eq(dbl, ggsw_dft->mat[(i - 1) * (K_TILDEBASE + 1) * NLIMBSBASE * NBASE + j * NLIMBSBASE * NBASE], 0.1));
			cr_assert(eq(dbl,
			             ggsw_dft->mat[(i - 1) * (K_TILDEBASE + 1) * NLIMBSBASE * NBASE + j * NLIMBSBASE * NBASE + 1],
			             0.2));
		}

	// Clean up
	delete_ggsw_dft(ggsw_dft);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
}

Test(add_ggsw_dft, basic)
{
	// Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	GGSWCiphertextDFT* ggsw_lhs_dft     = new_ggsw_dft(params_ggsw);
	GGSWCiphertextDFT* ggsw_rhs_dft     = new_ggsw_dft(params_ggsw);
	GGSWCiphertextDFT* sum_computed_dft = new_ggsw_dft(params_ggsw);

	// Draws uniformly the bivGGSW ciphertexts
	uniform_random_vec_znx_dft(module, ggsw_lhs_dft->mat, ggsw_size(params_ggsw), KAPPABASE - 1);
	uniform_random_vec_znx_dft(module, ggsw_rhs_dft->mat, ggsw_size(params_ggsw), KAPPABASE - 1);

	// Computes ggsw_lhs_dft + ggsw_rhs_dft
	add_ggsw_dft(sum_computed_dft, ggsw_lhs_dft, ggsw_rhs_dft);

	uint64_t nb_rows = params_ggsw->n_limbs_tilde;
	uint64_t nb_cols = params_glwe->n_limbs;
	uint64_t N       = params_glwe->nn;

	// Asserts sum_computed = ggsw_lhs_dft + ggsw_rhs_dft
	for (uint64_t idx = 0; idx < N * nb_cols * nb_rows; ++idx)
		cr_assert(eq(dbl, sum_computed_dft->mat[idx], ggsw_lhs_dft->mat[idx] + ggsw_rhs_dft->mat[idx]),
		          "add_biv_ggswy mismatch at index %" PRId64 ": %" PRId64 " + %" PRId64 " = %" PRId64 ", got %" PRId64,
		          (long long)idx, ggsw_lhs_dft->mat[idx], ggsw_rhs_dft->mat[idx],
		          ggsw_lhs_dft->mat[idx] + ggsw_rhs_dft->mat[idx], sum_computed_dft->mat[idx]);

	// Clean up
	delete_ggsw_dft(ggsw_lhs_dft);
	delete_ggsw_dft(ggsw_rhs_dft);
	delete_ggsw_dft(sum_computed_dft);
	pvda_delete_module_info(module);
	delete_ggsw_params(params_ggsw);
	delete_glwe_params(params_glwe);
}

Test(const_mult_ggsw_dft, without_normalization)
{
	// Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	PolyUnivDFT* u_dft                   = new_univ_dft(module);
	GGSWCiphertextDFT* ggsw_dft          = new_ggsw_dft(params_ggsw);
	GGSWCiphertextDFT* prod_computed_dft = new_ggsw_dft(params_ggsw);
	PolyUniv* u                          = new_univ(params_glwe);
	GGSWCiphertext* ggsw_ct              = new_ggsw(params_ggsw);
	GGSWCiphertext* prod_comp            = new_ggsw(params_ggsw);
	PolyUniv* prod_expected              = new_univ(params_glwe);

	// Draws uniformly the Zn[X] polynomial in the DFT domain
	uniform_random_vec_znx_dft(module, u_dft, 1, KAPPABASE - 1);

	// Draws uniformly the bivGGSW ciphertext in the DFT domain
	uniform_random_vec_znx_dft(module, ggsw_dft->mat, ggsw_size(params_ggsw), KAPPABASE - 1);

	// Computes DFT(u) * DFT(ggsw)
	const_mult_ggsw_dft(module, prod_computed_dft, ggsw_dft, u_dft);

	// Computes the matrix of u_dft, ggsw_dft and prod_computed_dft out of the DFT domain
	pvda_vec_znx_idft(module, u, 1, u_dft, 1);
	pvda_vec_znx_idft(module, ggsw_ct->mat, ggsw_size(params_ggsw), ggsw_dft->mat, ggsw_size(params_ggsw));
	pvda_vec_znx_idft(module, prod_comp->mat, ggsw_size(params_ggsw), prod_computed_dft->mat, ggsw_size(params_ggsw));

	for (uint64_t ii = 1; ii <= ggsw_num_pggsw(params_ggsw); ii++)
		for (uint64_t jj = 0; jj < K_TILDEBASE + 1; jj++)
		{
			VecBiv* ct_mat_ii_jj  = ggsw_retrieve_bivglwe(ggsw_ct, jj, ii);
			VecBiv* res_mat_ii_jj = ggsw_retrieve_bivglwe(prod_comp, jj, ii);
			for (uint64_t j = 0; j < KBASE + 1; j++)
				for (uint64_t i = LBASE; i >= 1; i--)
				{
					PolyUniv* ct_ij = ct_mat_ii_jj + (i - 1) * (KBASE + 1) * NBASE + j * NBASE;
					pvda_znx_product(module, prod_expected, u, ct_ij);
					for (uint64_t p = 0; p < NBASE; p++)
					{
						cr_assert(
						    eq(i64, res_mat_ii_jj[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p], prod_expected[p]));
					}
				}
		}

	// Clean up
	delete_univ(u);
	delete_univ_dft(u_dft);
	delete_univ(prod_expected);
	delete_ggsw(ggsw_ct);
	delete_ggsw(prod_comp);
	delete_ggsw_dft(ggsw_dft);
	delete_ggsw_dft(prod_computed_dft);
	delete_ggsw_params(params_ggsw);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}
