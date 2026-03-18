#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "ggsw_params.h"

#define NBASE            1024
#define KBASE            8
#define KAPPABASE        4
#define NLIMBSBASE       45
#define LBASE            NLIMBSBASE / (KBASE + 1)
#define SIGMABASE        1e-7

#define K_TILDEBASE      1
#define KAPPA_TILDEBASE  4
#define NLIMBS_TILDEBASE 10
#define L_TILDEBASE      NLIMBS_TILDEBASE / (K_TILDEBASE + 1)

//! COMMON PART (begin)

//! bivGGSW Part (begin)

// Test normalize_ggsw
Test(new_ggsw_params, basic)
{
	GLWEParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	cr_assert(params_ggsw != NULL);

	delete_glwe_ct_params(params_glwe);
	delete_ggsw_ct_params(params_ggsw);
}
