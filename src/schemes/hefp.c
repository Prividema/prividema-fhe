#include "hefp.h"

#include <assert.h>
#include <stdint.h>

#include "backend.h"
#include "backend_private.h"  //todo remove
#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "hefp_private.h"
#include "univariate_polynomial.h"
#include "utils.h"

int hefp_encode(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* out, uint64_t n_slots,
                int64_t scale_bits, double complex* in)
{
	int status = -1;

	PolyUnivRnX* tmp = new_univ_rnx(params_glwe);

	CHECK_ALLOC(tmp, "Tmp vector alloc failed in HE Fixed Point encoding");

	CHECK_CALL(hefp_encode_internal(module, tmp, n_slots, in, 0), "HE Fixed Point encoding failed");

	CHECK_CALL(univ_rnx_to_biv(params_glwe, out, tmp, -scale_bits),
	           "Biv conversion + scaling in HE Fixed Point encoding failed");

	status = 0;
cleanup:

	delete_univ_rnx(tmp);
	return status;
}

int hefp_decode(const PvdaBackend* module, const GLWEParams* params_glwe, double complex* out, uint64_t n_slots,
                int64_t scale_bits, PolyBiv* in)

{
	int status = -1;

	PolyUnivRnX* tmp = new_univ_rnx(params_glwe);

	CHECK_ALLOC(tmp, "Tmp vector alloc failed in HE Fixed Point decoding");

	biv_to_univ_rnx(params_glwe, tmp, in, -scale_bits);

	CHECK_CALL(hefp_decode_internal(module, out, n_slots, tmp), "HE Fixed Point decoding failed");

	status = 0;
cleanup:

	delete_univ_rnx(tmp);
	return status;
}
