#ifndef PVDA_FFTS_H
#define PVDA_FFTS_H

#include <complex.h>
#include <stdint.h>

struct pvda_fft_data_t
{
	uint64_t nn;
	double _Complex* roots;
	uint32_t* rotation_group;
};

struct pvda_fft_data_t* new_fft_data(uint64_t nn);

void delete_fft_data(struct pvda_fft_data_t* data);

#endif
