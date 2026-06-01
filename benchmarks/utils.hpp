#ifndef PVDA_BM_UTILS
#define PVDA_BM_UTILS

extern "C" {
#include "glwe_params.h"
#include "univariate_polynomial.h"
};

int rnx_random_vec(PolyUnivRnX* res, GLWEParams* params_glwe);

#endif  //PVDA_BM_UTILS
