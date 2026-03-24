#include "maths_structures.h"

#include <stdio.h>

void printf_poly_biv(PolyBiv* pol, int64_t pol_sl, int64_t nn, int64_t l)
{
	printf("\n");
	for (uint64_t i = 1; i <= l; i++)
	{
		printf("\nY^%ld : ", i);
		for (uint64_t p = 0; p < nn; p++)
			printf((pol[(i - 1) * pol_sl + p] < 0) ? "%ld X^%ld " : "%ld X^%ld ", pol[(i - 1) * pol_sl + p], p);
	}
}

void printf_vec_poly_biv(VecBiv* pols, int64_t pols_nb, int64_t nn, int64_t l)
{
	printf("\n");
	for (uint64_t j = 0; j < pols_nb; j++)
	{
		printf("\n%ld-th component : ", j);
		printf_poly_biv(pols + j * nn, pols_nb * nn, nn, l);
		printf("\n");
	}
}

void printf_poly_univ_ZnX(PolyUniv* pol, int64_t nn)
{
	printf((pol[0] < 0) ? "%ld " : " %ld ", pol[0]);
	for (uint64_t p = 1; p < nn; p++)
		printf((pol[p] < 0) ? "- %ld X^%ld " : "+ %ld X^%ld ", (pol[p] < 0) ? -pol[p] : pol[p], p);
}

void printf_poly_univ_RnX(double* pol, int64_t nn)
{
	for (uint64_t p = 0; p < nn; p++) printf(" %lf X^%ld ", pol[p], p);
}

void printf_vec_poly_univ(VecBiv* pols, int64_t pols_size, int64_t nn)
{
	printf("\n");
	for (uint64_t j = 0; j < pols_size; j++)
	{
		printf("\n%ld-th component : ", j);
		printf_poly_univ_ZnX(pols + j * nn, nn);
		printf("\n");
	}
}

void printf_secret_key(PolyUniv** sk_values, int64_t nn, int64_t k)
{
	printf("\n\nBegin Secret Key:\n");
	for (uint64_t j = 0; j < k; j++)
	{
		printf("\n%ld-th component", j);
		printf_poly_univ_ZnX(sk_values[j], nn);
	}
	printf("\n\nEnd Secret Key\n");
}
