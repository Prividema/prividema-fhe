#ifndef PVDA_BM_UTILS
#define PVDA_BM_UTILS

extern "C" {
#include "backend.h"
#include "glwe_params.h"
#include "univariate_polynomial.h"
};

int rnx_random_vec(const PvdaBackend* module, PolyUnivRnX* res, GLWEParams* params_glwe);

#endif  //PVDA_BM_UTILS
