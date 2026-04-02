
#include "ututils.h"

#include <math.h>

#include "math.h"

double generate_sigma(PvdaTstParams* p)
{
	if (p->sigma > 0.0) return p->sigma;
	if (p->sigma < 0.0) return ldexp(1.0, (int)p->sigma);
	return ldexp(1.0, -(p->l / 2 + 1) * p->kappa);
}

/*

OLD CODE that might be useful later on

void printf_glwe(GLWECiphertext* glwe)
{
    for (int64_t j = 0; j < KBASE + 1; j++) printf_poly_biv(glwe->vec + j * NBASE, (KBASE + 1) * NBASE, NBASE, LBASE);
}

void printf_glwe_dft(MODULE* module, GLWECiphertextDFT* glwe_dft)
{
    VecBiv* glwe_vec = malloc(glwe_params_bytes(glwe_dft->params));
    pvda_vec_znx_idft(module, glwe_vec, glwe_params_n_limbs(glwe_dft->params), glwe_dft->vec,
                      glwe_params_n_limbs(glwe_dft->params));

    for (int64_t j = 0; j < KBASE + 1; j++) printf_poly_biv(glwe_vec + j * NBASE, (KBASE + 1) * NBASE, NBASE, LBASE);

    free(glwe_vec);
}

*/
