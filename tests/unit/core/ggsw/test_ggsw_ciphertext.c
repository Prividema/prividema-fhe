#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdint.h>
#include <sys/types.h>

#include "bivariate_polynomial.h"
#include "ggsw.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe_params.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

#define NLIMBSBASE       (params_glwe->k * glwe_params_l_a(params_glwe) + glwe_params_l_b(params_glwe))
#define NLIMBS_TILDEBASE (params_glwe->k * ggsw_params_l_tilde_a(params_ggsw) + ggsw_params_l_tilde_b(params_ggsw))

//! COMMON PART (begin)

// Test ggsw_size
PvdaParamTest(ggsw_size, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Asserts ggsw_size returns NLIMBS_TILDEBASE * NLIMBSBASE
	cr_assert(eq(i64, ggsw_total_n_glwe_limbs(params_ggsw), NLIMBS_TILDEBASE * NLIMBSBASE));

	DELETE_PVDA_PARAMS_GGSW;
}

// Test ggsw_bytes
PvdaParamTest(ggsw_bytes, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Asserts ggsw_bytes returns NLIMBS_TILDEBASE * NLIMBSBASE * params_glwe->nn * sizeof(int64_t)
	cr_assert(eq(i64, ggsw_bytes(params_ggsw), NLIMBS_TILDEBASE * NLIMBSBASE * params_glwe->nn * sizeof(int64_t)));

	DELETE_PVDA_PARAMS_GGSW;
}

//! bivGGSW Part (begin)

// Test ggsw_coef_number
PvdaParamTest(ggsw_coef_number, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Asserts ggsw_coeff_number returns NLIMBS_TILDEBASE * NLIMBSBASE * params_glwe->nn
	cr_assert(eq(i64, ggsw_coef_number(params_ggsw), NLIMBS_TILDEBASE * NLIMBSBASE * params_glwe->nn));

	DELETE_PVDA_PARAMS_GGSW;
}

// Test new_ggsw
PvdaParamTest(new_ggsw, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Variables
	GGSWCiphertext* ggsw = new_ggsw(params_ggsw);

	// Asserts new_ggsw worked
	cr_assert(eq(int, ggsw->mat != NULL, 1));
	cr_assert(eq(int, ggsw->params != NULL, 1));

	// Clean up
	delete_ggsw(ggsw);

	DELETE_PVDA_PARAMS_GGSW;
}

// Test ggsw_Sj_Yti
PvdaParamTest(ggsw_Sj_Yti, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Variables
	GGSWCiphertext* ggsw = new_ggsw(params_ggsw);

	// asserts ggsw_Sj_Yti returns the right pointer
	for (uint64_t ij = 0; ij < ggsw_num_rows(params_ggsw); ++ij)
	{
		uint64_t j        = (ij % (params_ggsw->k_tilde + 1));
		uint64_t i        = ij / (params_ggsw->k_tilde + 1) + 1;
		VecBiv* ct_mat_ij = ggsw_retrieve_bivglwe(ggsw, j, i);

		// Modify the two firsts coefficients of biGLWE(-m * sk_j / 2^{kappa_tilde * i}).
		ct_mat_ij[0] = 1;
		ct_mat_ij[1] = 2;

		cr_assert(eq(
		    i64,
		    ggsw->mat[(i - 1) * (params_glwe->k + 1) * NLIMBSBASE * params_glwe->nn + j * NLIMBSBASE * params_glwe->nn],
		    1));
		cr_assert(eq(i64,
		             ggsw->mat[(i - 1) * (params_glwe->k + 1) * NLIMBSBASE * params_glwe->nn +
		                       j * NLIMBSBASE * params_glwe->nn + 1],
		             2));
	}

	// Clean up
	delete_ggsw(ggsw);

	DELETE_PVDA_PARAMS_GGSW;
}

// Test normalize_ggsw
PvdaParamTest(normalize_ggsw, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	GGSWCiphertext* ggsw = new_ggsw(params_ggsw);
	GGSWCiphertext* res  = new_ggsw(params_ggsw);

	// Normalize the ggsw ciphertext
	int status = normalize_ggsw(module, res, ggsw);

	// Asserts normalize_ggsw succeed
	cr_assert(eq(int, status, 0), "normalize_ggsw failed");

	// Clean up
	delete_ggsw(ggsw);
	delete_ggsw(res);

	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(add_ggsw, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Variables
	GGSWCiphertext* ggsw_lhs     = new_ggsw(params_ggsw);
	GGSWCiphertext* ggsw_rhs     = new_ggsw(params_ggsw);
	GGSWCiphertext* sum_computed = new_ggsw(params_ggsw);

	// Draws uniformly the bivGGSW ciphertexts
	uniform_random_vec(params_glwe->nn, ggsw_lhs->mat, ggsw_total_n_glwe_limbs(params_ggsw), params_glwe->nn,
	                   params_glwe->kappa - 1);
	uniform_random_vec(params_glwe->nn, ggsw_rhs->mat, ggsw_total_n_glwe_limbs(params_ggsw), params_glwe->nn,
	                   params_glwe->kappa - 1);

	// Computes ggsw_lhs + ggsw_rhs
	add_ggsw(sum_computed, ggsw_lhs, ggsw_rhs);

	uint64_t nb_rows = ggsw_num_rows(sum_computed->params);
	uint64_t nb_cols = glwe_params_n_limbs(sum_computed->params->params_glwe);
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

	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(const_mult_ggsw, without_normalization, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Variables
	PolyUniv* u                      = new_univ(params_glwe);
	PolyUnivDFT* u_dft               = new_univ_dft(module);
	GGSWCiphertext* ggsw             = new_ggsw(params_ggsw);
	GGSWCiphertext* product_computed = new_ggsw(params_ggsw);
	PolyUniv* prod_expected          = new_univ(params_glwe);

	// Draws uniformly the Zn[X] constant
	uniform_random_vec(params_glwe->nn, u, 1, params_glwe->nn, params_glwe->kappa - 1);

	// Draws uniformly the bivGGSW ciphertext
	uniform_random_vec(params_glwe->nn, ggsw->mat, ggsw_total_n_glwe_limbs(params_ggsw), params_glwe->nn,
	                   params_glwe->kappa - 1);

	// Computes u in the DFT domain
	univ_coefs_to_dft(module, u_dft, u);

	// Computes u * ggsw
	const_mult_ggsw(module, product_computed, ggsw, u_dft);

	// Asserts product_computed = u * ggsw
	for (uint64_t iijj = 0; iijj < ggsw_num_rows(params_ggsw); ++iijj)
	{
		uint64_t jj           = (iijj % (params_ggsw->k_tilde + 1));
		uint64_t ii           = iijj / (params_ggsw->k_tilde + 1) + 1;
		VecBiv* ct_mat_ii_jj  = ggsw_retrieve_bivglwe(ggsw, jj, ii);
		VecBiv* res_mat_ii_jj = ggsw_retrieve_bivglwe(product_computed, jj, ii);
		for (uint64_t ij = 0; ij < glwe_params_n_limbs(params_glwe); ++ij)
		{
			uint64_t j      = (ij % (params_glwe->k + 1));
			uint64_t i      = ij / (params_glwe->k + 1) + 1;
			PolyUniv* ct_ij = ct_mat_ii_jj + (i - 1) * (params_glwe->k + 1) * params_glwe->nn + j * params_glwe->nn;
			pvda_znx_small_product(module, prod_expected, u, ct_ij);
			for (uint64_t p = 0; p < params_glwe->nn; p++)
			{
				cr_assert(eq(i64,
				             res_mat_ii_jj[(i - 1) * (params_glwe->k + 1) * params_glwe->nn + j * params_glwe->nn + p],
				             prod_expected[p]));
			}
		}
	}

	// Clean up
	delete_univ(u);
	delete_univ(prod_expected);
	delete_univ_dft(u_dft);
	delete_ggsw(ggsw);
	delete_ggsw(product_computed);

	DELETE_PVDA_PARAMS_GGSW;
}

//! bivGGSW IN DFT SPACE Part (begin)

// Test ggsw_coef_number
PvdaParamTest(ggsw_coef_number_dft, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Asserts ggsw_coef_number_dft returns NLIMBS_TILDEBASE * NLIMBSBASE * params_glwe->nn / 2
	cr_assert(eq(i64, ggsw_coef_number_dft(params_ggsw), NLIMBS_TILDEBASE * NLIMBSBASE * params_glwe->nn / 2));

	DELETE_PVDA_PARAMS_GGSW;
}

// Test new_ggsw
PvdaParamTest(new_ggsw_dft, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Variables
	GGSWCiphertextDFT* ggsw_dft = new_ggsw_dft(params_ggsw);

	// Asserts new_ggsw_dft worked
	cr_assert(eq(int, ggsw_dft->mat != NULL, 1));
	cr_assert(eq(int, ggsw_dft->params != NULL, 1));

	// Clean up
	delete_ggsw_dft(ggsw_dft);

	DELETE_PVDA_PARAMS_GGSW;
}

// Test ggsw_Sj_Yti_dft
PvdaParamTest(ggsw_Sj_Yti_dft, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Variables
	GGSWCiphertextDFT* ggsw_dft = new_ggsw_dft(params_ggsw);

	// Asserts ggsw_Sj_Yti_dft returns the right pointer
	for (uint64_t ij = 0; ij < ggsw_num_rows(params_ggsw); ++ij)
	{
		uint64_t j           = (ij % (params_ggsw->k_tilde + 1));
		uint64_t i           = ij / (params_ggsw->k_tilde + 1) + 1;
		VecBivDFT* ct_mat_ij = ggsw_retrieve_bivglwe_dft(ggsw_dft, j, i);

		// Modify the two firsts coefficients of biGLWE(-m * sk_j / 2^{kappa_tilde * i}).
		ct_mat_ij[0] = 0.1;
		ct_mat_ij[1] = 0.2;

		cr_assert(eq(dbl,
		             ggsw_dft->mat[(i - 1) * (params_glwe->k + 1) * NLIMBSBASE * params_glwe->nn +
		                           j * NLIMBSBASE * params_glwe->nn],
		             0.1));
		cr_assert(eq(dbl,
		             ggsw_dft->mat[(i - 1) * (params_glwe->k + 1) * NLIMBSBASE * params_glwe->nn +
		                           j * NLIMBSBASE * params_glwe->nn + 1],
		             0.2));
	}

	// Clean up
	delete_ggsw_dft(ggsw_dft);

	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(add_ggsw_dft, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Variables
	GGSWCiphertextDFT* ggsw_lhs_dft     = new_ggsw_dft(params_ggsw);
	GGSWCiphertextDFT* ggsw_rhs_dft     = new_ggsw_dft(params_ggsw);
	GGSWCiphertextDFT* sum_computed_dft = new_ggsw_dft(params_ggsw);

	// Draws uniformly the bivGGSW ciphertexts
	uniform_random_vec_znx_dft(module, ggsw_lhs_dft->mat, ggsw_total_n_glwe_limbs(params_ggsw), params_glwe->kappa - 1);
	uniform_random_vec_znx_dft(module, ggsw_rhs_dft->mat, ggsw_total_n_glwe_limbs(params_ggsw), params_glwe->kappa - 1);

	// Computes ggsw_lhs_dft + ggsw_rhs_dft
	add_ggsw_dft(sum_computed_dft, ggsw_lhs_dft, ggsw_rhs_dft);

	uint64_t nb_rows = ggsw_num_rows(params_ggsw);
	uint64_t nb_cols = glwe_params_n_limbs(params_glwe);
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

	DELETE_PVDA_PARAMS_GGSW;
}

PvdaParamTest(const_mult_ggsw_dft, without_normalization, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	// Variables
	PolyUnivDFT* u_dft                   = new_univ_dft(module);
	GGSWCiphertextDFT* ggsw_dft          = new_ggsw_dft(params_ggsw);
	GGSWCiphertextDFT* prod_computed_dft = new_ggsw_dft(params_ggsw);
	PolyUniv* u                          = new_univ(params_glwe);
	GGSWCiphertext* ggsw_ct              = new_ggsw(params_ggsw);
	GGSWCiphertext* prod_comp            = new_ggsw(params_ggsw);
	PolyUniv* prod_expected              = new_univ(params_glwe);

	// Draws uniformly the Zn[X] polynomial in the DFT domain
	uniform_random_vec_znx_dft(module, u_dft, 1, params_glwe->kappa - 1);

	// Draws uniformly the bivGGSW ciphertext in the DFT domain
	uniform_random_vec_znx_dft(module, ggsw_dft->mat, ggsw_total_n_glwe_limbs(params_ggsw), params_glwe->kappa - 1);

	// Computes DFT(u) * DFT(ggsw)
	const_mult_ggsw_dft(module, prod_computed_dft, ggsw_dft, u_dft);

	// Computes the matrix of u_dft, ggsw_dft and prod_computed_dft out of the DFT domain
	pvda_vec_znx_idft(module, u, 1, u_dft, 1);
	pvda_vec_znx_idft(module, ggsw_ct->mat, ggsw_total_n_glwe_limbs(params_ggsw), ggsw_dft->mat,
	                  ggsw_total_n_glwe_limbs(params_ggsw));
	pvda_vec_znx_idft(module, prod_comp->mat, ggsw_total_n_glwe_limbs(params_ggsw), prod_computed_dft->mat,
	                  ggsw_total_n_glwe_limbs(params_ggsw));

	for (uint64_t iijj = 0; iijj < ggsw_num_rows(params_ggsw); ++iijj)
	{
		uint64_t jj           = (iijj % (params_ggsw->k_tilde + 1));
		uint64_t ii           = iijj / (params_ggsw->k_tilde + 1) + 1;
		VecBiv* ct_mat_ii_jj  = ggsw_retrieve_bivglwe(ggsw_ct, jj, ii);
		VecBiv* res_mat_ii_jj = ggsw_retrieve_bivglwe(prod_comp, jj, ii);
		for (uint64_t ij = 0; ij < glwe_params_n_limbs(params_glwe); ++ij)
		{
			uint64_t j      = ij % (params_glwe->k + 1);
			uint64_t i      = ij / (params_glwe->k + 1) + 1;
			PolyUniv* ct_ij = ct_mat_ii_jj + (i - 1) * (params_glwe->k + 1) * params_glwe->nn + j * params_glwe->nn;
			pvda_znx_small_product(module, prod_expected, u, ct_ij);
			for (uint64_t p = 0; p < params_glwe->nn; p++)
			{
				cr_assert(eq(i64,
				             res_mat_ii_jj[(i - 1) * (params_glwe->k + 1) * params_glwe->nn + j * params_glwe->nn + p],
				             prod_expected[p]));
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

	DELETE_PVDA_PARAMS_GGSW;
}
