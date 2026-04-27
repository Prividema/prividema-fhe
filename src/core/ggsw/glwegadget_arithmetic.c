#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "glwegadget_arithmetic.h"
#include "glwegadget_ciphertext.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

int glwegadget_half_prod(const MODULE* module, GLWECiphertext* result,
                         const GLWEGadgetCiphertextPrep* glwegadget_prep_ct, const PolyBiv* a)
{
	int status = -1;

	//TODO: assert size compatibility
	size_t nrows = glwegadget_prep_ct->params->l_tilde;
	uint64_t nn  = glwegadget_prep_ct->params->params_glwe->nn;
	size_t ncols = glwe_params_n_limbs(glwegadget_prep_ct->params->params_glwe);

	GLWECiphertextDFT* glwe_dft = new_glwe_dft(result->params);

	CHECK_CALL(pvda_vmp_apply_dft(module, glwe_dft->vec, ncols, a, nrows, nn, glwegadget_prep_ct->mat, nrows, ncols),
	           "vmp apply falied in half product");

	CHECK_CALL(glwe_dft_to_coef(module, result, glwe_dft),
	           "conversion from GLWE DFT to coefs failed in GLWEGadget half product");
	status = 0;
cleanup:
	delete_glwe_dft(glwe_dft);

	return status;
}

int glwegadget_half_prod_dft_to_dft(const MODULE* module, GLWECiphertextDFT* result_dft,
                                    const GLWEGadgetCiphertextPrep* glwegadget_prep_ct, const PolyBivDFT* a_dft)
{
	int status = -1;

	//TODO: assert size compatibility
	size_t nrows = glwegadget_prep_ct->params->l_tilde;
	size_t ncols = glwe_params_n_limbs(glwegadget_prep_ct->params->params_glwe);

	CHECK_CALL(pvda_vmp_apply_dft_to_dft(module, result_dft->vec, ncols, (double*)a_dft, nrows, glwegadget_prep_ct->mat,
	                                     nrows, ncols),
	           "vmp apply falied in half product");

	status = 0;
cleanup:
	return status;
}
