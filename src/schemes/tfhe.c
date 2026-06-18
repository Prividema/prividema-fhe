#include "tfhe.h"

#include <string.h>

#include "ggsw_arithmetic.h"
#include "glwe_arithmetic.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "utils.h"

int tfhe_cmux_unprepared(const MODULE* module, GLWECiphertext* res, const GLWECiphertext* c0, const GLWECiphertext* c1,
                         const GGSWCiphertext* sel, int normalize_sub)
{
	//TODO: this only supports res of equal precision than the inputs at the moment
	// Since there is an external product, maybe at some point we want to drop some precision?
	sub_glwe(module, res, c1, c0);
	if (normalize_sub) CHECK_CALL(normalize_glwe(module, res, res), "normalization failed in CMux");
	CHECK_CALL(ggsw_unprepared_external_product(module, res, res, sel), "GGSW external product failed in CMux");
	add_glwe(module, res, res, c0);

	return 0;
cleanup:
	return -1;
}

int tfhe_cmux(const MODULE* module, GLWECiphertext* res, const GLWECiphertext* c0, const GLWECiphertext* c1,
              const GGSWCiphertextPrep* sel, int normalize_sub)
{
	//TODO: this only supports res of equal precision than the inputs at the moment
	// we might want to allow a precision drop
	sub_glwe(module, res, c1, c0);
	if (normalize_sub) CHECK_CALL(normalize_glwe(module, res, res), "normalization failed in CMux");
	CHECK_CALL(ggsw_external_product(module, res, res, sel), "GGSW external product failed in CMux");
	add_glwe(module, res, res, c0);

	return 0;
cleanup:
	return -1;
}

int tfhe_cmux_tree(const MODULE* module, GLWECiphertext* res, const GLWECiphertext** src, int inp_cols,
                   const GGSWCiphertextPrep** selectors, int sel_size, int delete_src)
{
	GLWEParams* aggregation_params = res->params;
	int64_t log_inp_cols           = next_pow2_log(inp_cols);

	GLWECiphertext* glwe_tree[(log_inp_cols + 1)][inp_cols];  //Inefficient, suffices for now
	memset(glwe_tree, 0, (log_inp_cols + 1) * inp_cols * sizeof(GLWECiphertext*));

	//Note: this check MUST go after glwe_tree definition and set-to-0 for correct error handling
	if (sel_size < log_inp_cols)
	{
		RAISE_ERROR("Tried to use a CMux tree with less selection signals than the required for the number of inputs");
	}

	for (int c = 0; c < inp_cols; ++c)
	{
		glwe_tree[0][c] = src[c];
	}

	int used_cols = inp_cols;
	for (int l = 0; l < log_inp_cols; ++l)
	{
		//Init next level
		if (l == log_inp_cols - 1)
			glwe_tree[l + 1][0] = res;
		else
			for (int c = 0; c < (used_cols + 1) / 2; ++c)
			{
				glwe_tree[l + 1][c] = new_glwe(aggregation_params);
				CHECK_ALLOC(glwe_tree[l + 1][c], "GLWE allocation failed in a CMux tree");
			}

		//CMux
		int c;
		for (c = 0; c < used_cols; c += 2)
		{
			if (c + 1 < used_cols)
				CHECK_CALL(
				    tfhe_cmux(module, glwe_tree[l + 1][c / 2], glwe_tree[l][c], glwe_tree[l][c + 1], selectors[l], 1),
				    "CMux failed in a CMux tree");
			else
			{
				glwe_copy(glwe_tree[l + 1][c / 2], glwe_tree[l][c]);
			}
		}

		for (int c = 0; c < used_cols; ++c)
		{
			if (l != 0 || delete_src)
			{
				delete_glwe(glwe_tree[l][c]);
				glwe_tree[l][c] = 0;
			}
		}
		used_cols = (used_cols + 1) / 2;
	}

	return 0;
cleanup:
	//Cleanup logic depends on whether we want to remove the first row or not
	for (int l = delete_src ? 0 : 1; l < log_inp_cols; ++l)
	{
		for (int c = 0; c < inp_cols; ++c) delete_glwe(glwe_tree[l][c]);
	}
	return -1;
}
