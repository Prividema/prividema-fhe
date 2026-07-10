#include "pvda_ffts.h"

#include <math.h>
#include <stdlib.h>

#include "utils.h"

static inline void init_roots(double complex* vec, int nn)
{
	double complex r = cexp(M_PI * I / nn);
	double complex v = 1;
	for (int i = 0; i < 2 * nn + 1; ++i)
	{
		vec[i] = v;
		v *= r;
	}
}

static inline void init_group(uint32_t* vec, int nn)
{
	uint32_t r    = 5;
	uint32_t v    = 1;
	uint32_t mask = 2 * nn - 1;

	for (int i = 0; i < nn / 4; ++i)
	{
		vec[i] = v;
		v *= r;
		v &= mask;
	}
}

struct pvda_fft_data_t* new_fft_data(uint64_t nn)
{
	struct pvda_fft_data_t* ans = malloc(sizeof(struct pvda_fft_data_t));
	CHECK_ALLOC(ans, "failed allocation of Prividema FFT data");

	ans->nn             = nn;
	ans->roots          = calloc(2 * nn + 1, sizeof(double complex));
	ans->rotation_group = calloc(nn / 4, sizeof(int32_t));
	CHECK_ALLOC(ans->roots, "failed allocation of roots in Privideam FFT data");
	CHECK_ALLOC(ans->rotation_group, "failed allocation of rotation group in Privideam FFT data");

	init_roots(ans->roots, nn);
	init_group(ans->rotation_group, nn);

	return ans;
cleanup:
	delete_fft_data(ans);
	return NULL;
}

void delete_fft_data(struct pvda_fft_data_t* data)
{
	if (!data) return;
	free(data->roots);
	free(data->rotation_group);
	free(data);
}
