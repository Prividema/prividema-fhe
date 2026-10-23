#include "pvda_ffts.h"

#include <math.h>
#include <stdlib.h>

#include "backend.h"
#include "backend_nt.h"
#include "backend_private.h"
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

	for (int i = 0; i < nn; ++i)
	{
		vec[i] = v;
		v *= r;
		v &= mask;
	}
}

struct pvda_fft_data_t* new_fft_data(uint64_t nn)
{
	// Maximum number of NTT tables that the struct will be able to hold
	const int num_ntt_tables = 10;

	struct pvda_fft_data_t* ans = malloc(sizeof(struct pvda_fft_data_t));
	CHECK_ALLOC(ans, "failed allocation of Prividema FFT data");

	ans->nn             = nn;
	ans->roots          = calloc(2 * nn + 1, sizeof(double complex));
	ans->rotation_group = calloc(nn, sizeof(int32_t));
	ans->num_ntt_tables = num_ntt_tables;
	ans->ntt_tables     = calloc(num_ntt_tables, sizeof(struct ntt_root_table_t));
	CHECK_ALLOC(ans->roots, "failed allocation of roots in Prividema FFT data");
	CHECK_ALLOC(ans->rotation_group, "failed allocation of rotation group in Prividema FFT data");
	CHECK_ALLOC(ans->ntt_tables, "failed allocation for NTT tables in prividema FFT data");

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

	if (data->ntt_tables)
		for (uint64_t i = 0; i < data->num_ntt_tables; ++i)
		{
			free(data->ntt_tables[i].roots);
		}
	free(data->ntt_tables);

	free(data);
}

NTTRoot* fft_data_get_ntt_table(struct pvda_fft_data_t* data, uint64_t t)
{
	for (uint64_t i = 0; i < data->num_ntt_tables; ++i)
	{
		if (data->ntt_tables[i].t == t) return data->ntt_tables[i].roots;
	}
	return NULL;
}

int ntt_table_generate(const PvdaBackend* module, uint64_t t)
{
	struct ntt_root_table_t* table_entry = NULL;
	struct pvda_fft_data_t* data         = module->pvda_fft_data;
	for (uint64_t i = 0; i < data->num_ntt_tables; ++i)
	{
		if (data->ntt_tables[i].t == t)
		{
			return 1;
		}
		if (data->ntt_tables[i].t == 0)
		{
			table_entry = &data->ntt_tables[i];
		}
	}

	if (!table_entry) RAISE_ERROR("No free NTT table slot left in Prividema FFT data");

	const uint64_t n_elems = 2 * data->nn + 1;

	table_entry->roots = calloc(n_elems, sizeof(NTTRoot));
	CHECK_ALLOC(table_entry->roots, "failed allocation of NTT root table");
	table_entry->t = t;

	NTTRoot r = pvda_ring_2nth_root(module, t);
	NTTRoot v = 1;

	for (uint64_t i = 0; i < n_elems; ++i)
	{
		table_entry->roots[i] = v;
		v *= r;
		v %= t;
	}

	return 0;

cleanup:
	return -1;
}
