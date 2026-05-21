#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "univariate_polynomial.h"

void biv_to_univ_rnx_new(const GLWEParams* params, double val, PolyBiv* biv);
int main()
{
	printf("Hello example private match!\n");

	GLWEParams* params_glwe = new_glwe_params(1, 1, 19, 2 * 400, 0, NOISE_FAST_UNIFORM);

	PolyUnivRnX* rnx_final = new_univ_rnx(params_glwe);
	PolyBiv* biv           = new_biv_poly(params_glwe);

	/*
	  uint64_t stnx = (1UL << 62) + (1UL << 49) + 1 + 3000;
	  int exp       = 2;
	  double mag    = ldexp(stnx & ((1UL << 63) - 1), -exp);
	  if (mag > 1) mag -= (int)mag;
	  if (mag > 0.5) mag -= 1;
	  printf("expt: %.16e\n", mag);
	  printf("expt_T: %.16e\n", mag - (int)(mag));
	  _biv_decomp_internal(stnx, exp, biv, 1, params_glwe);
	*/

	double j = 4.125364e-322;
	biv_to_univ_rnx_new(params_glwe, j, biv);

	for (int i = 0; i < glwe_params_l_a(params_glwe); ++i)
	{
		printf("%ld ", biv[i]);
	}
	printf("\n");

	biv_to_univ_rnx(params_glwe, rnx_final, biv);
	double final = rnx_final[0];
	if (final > 0.5) final -= 1;
	printf("%.16e\n", final);

	delete_biv(biv);
	delete_univ_rnx(rnx_final);

	delete_glwe_params(params_glwe);

	return 0;
}
