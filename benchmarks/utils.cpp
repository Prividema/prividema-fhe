#include "utils.hpp"

#include <stdlib.h>

extern "C" {
#include "backend.h"
#include "rng.h"
};

int rnx_random_vec(const PvdaBackend* module, PolyUnivRnX* res, GLWEParams* params_glwe)

{
	int status = -1;

	PolyUnivTnX* tmp_tnx = new_univ_tnx(params_glwe);

	// (ab)use the fact that tnx and Z mod 2^64 are isomorphic and
	// the memory representation is the same for isomprphic values
	uniform_pow2_random_pol_znx(module, (PolyUniv*)tmp_tnx, params_glwe->nn, 64);

	univ_tnx_to_rnx(params_glwe, res, tmp_tnx);

	status = 0;
cleanup:
	delete_univ_tnx(tmp_tnx);
	return status;
}
