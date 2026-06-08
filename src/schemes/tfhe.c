#include "tfhe.h"

#include "ggsw_arithmetic.h"
#include "glwe_arithmetic.h"
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
	// Since there is an external product, maybe at some point we want to drop some precision?
	sub_glwe(module, res, c1, c0);
	if (normalize_sub) CHECK_CALL(normalize_glwe(module, res, res), "normalization failed in CMux");
	CHECK_CALL(ggsw_external_product(module, res, res, sel), "GGSW external product failed in CMux");
	add_glwe(module, res, res, c0);

	return 0;
cleanup:
	return -1;
}
