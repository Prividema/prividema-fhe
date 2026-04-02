#ifndef PVDA_UTEST_UTILS
#define PVDA_UTEST_UTILS

#include <stdint.h>
#include <sys/types.h>

int pvda_tst_util_mock();

typedef struct pvda_tst_params_t
{
	uint64_t nn;
	uint64_t k;
	uint64_t kappa;
	uint64_t l;
	uint64_t l_tilde;
	double sigma;
} PvdaTstParams;

double generate_sigma(PvdaTstParams* p);

#define INIT_PVDA_PARAMS_GLWE(PRS)                             \
	MODULE* module          = pvda_new_module_info((PRS)->nn); \
	double sigma            = generate_sigma((PRS));           \
	GLWEParams* params_glwe = new_glwe_params((PRS)->nn, (PRS)->k, (PRS)->kappa, (PRS)->l * ((PRS)->k + 1), sigma);

#define DELETE_PVDA_PARAMS_GLWE      \
	pvda_delete_module_info(module); \
	delete_glwe_params(params_glwe);

#define INIT_PVDA_PARAMS_GGSW(PRS) \
	INIT_PVDA_PARAMS_GLWE((PRS))   \
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, (PRS)->k, (PRS)->kappa, (PRS)->l_tilde * ((PRS)->k + 1))

#define DELETE_PVDA_PARAMS_GGSW \
	DELETE_PVDA_PARAMS_GLWE     \
	delete_ggsw_params(params_ggsw);

#endif
