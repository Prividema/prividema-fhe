#include "maths_structures.h"

#include <stdio.h>

void printf_poly_biv(PolyBiv* pol)
{
	printf("\n");
	for (uint64_t i = 0; i < pol->l; i++)
	{
		printf("\nY^%ld : ", i);
		for (uint64_t p = 0; p < pol->nn; p++) printf("%ld X^%ld ", pol->ptr[(i - 1) * pol->stride + p], p);
	}
}

void print_coefs_biv(const PolyBiv* biv, int max_n, int max_l)
{
	for (int l = 0; l < max_l && l < biv->l; ++l)
	{
		printf("Biv lvl %03d: ", l);

		for (int p = 0; p < max_n && p < biv->nn; ++p)
		{
			printf("%03ld ", biv->ptr[l * biv->stride + p]);
		}
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
