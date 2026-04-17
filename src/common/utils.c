#include "utils.h"

#include <math.h>
#include <stdint.h>

double torus_distance(double a, double b)
{
	a -= floor(a);
	b -= floor(b);
	if (a > b)
	{
		return fmin(a - b, (b + 1) - a);
	}
	return fmin(b - a, (a + 1) - b);
}

uint64_t u64_round_up_div(uint64_t num, uint64_t div) { return (num + (div - 1)) / div; }
int32_t i32_round_up_div(int32_t num, int32_t div) { return (num + (div - 1)) / div; }
