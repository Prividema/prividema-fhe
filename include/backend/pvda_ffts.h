#ifndef PVDA_FFTS_H
#define PVDA_FFTS_H

#include <complex.h>
#include <stdint.h>

typedef __int128_t NTTRoot;

struct ntt_root_table_t
{
	uint64_t t;
	NTTRoot* roots;
};

struct pvda_fft_data_t
{
	uint64_t nn;
	double _Complex* roots;
	uint32_t* rotation_group;

	uint64_t num_ntt_tables;
	struct ntt_root_table_t* ntt_tables;
};

struct pvda_fft_data_t* new_fft_data(uint64_t nn);

void delete_fft_data(struct pvda_fft_data_t* data);

/**
 *
 * Generates and inserts the NTT root table for the specified modulus t
 *
 * @retval 0 on success
 * @retval -1 on failure
 * @retval 1 if t is already present
 *
 */
int generate_ntt_table(struct pvda_fft_data_t* data, uint64_t t);

/**
 *
 *
 *
 *
 *
 */
NTTRoot* get_ntt_table(struct pvda_fft_data_t* data, uint64_t t);

#endif
