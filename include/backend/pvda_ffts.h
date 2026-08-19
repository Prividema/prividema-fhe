#ifndef PVDA_FFTS_H
#define PVDA_FFTS_H

#include <complex.h>
#include <stdint.h>

#include "backend.h"

typedef int64_t NTTRoot;

/**
 * @param Lists of 2N + 1 first powers of a primitive 2N-th root of unity
 *
 *
 */
struct ntt_root_table_t
{
	uint64_t t;      ///<  Modulus for which this list has been computed
	NTTRoot* roots;  ///< The list of powers of the primitive 2N-th root of unity modulo t
};

/**
 * @brief Struct for Prividema FFT data.
 *
 */
struct pvda_fft_data_t
{
	uint64_t nn;               ///< Parameter N
	double _Complex* roots;    ///< 2*N + 1 first powers of the 2Nth root of unity
	uint32_t* rotation_group;  ///<  List of powers of the form 5^j mod 2N,
	                           ///<  N first values

	uint64_t num_ntt_tables;              ///< length of ntt-tables
	struct ntt_root_table_t* ntt_tables;  ///< vector of root tables with associated modulus information
};

/**
 * @brief create a new Prividema-FHE FFT data object
 *
 * @param nn The parameter N (number of coefficients per polynomial)
 *
 * @return NULL in case of failure, a newly allocated object('s pointer') otherwise
 */
struct pvda_fft_data_t* new_fft_data(uint64_t nn);

/**
 * @brief Deletes Prividema-FHE FFT data
 *
 * @param data The Prividema-FHE FFT data object to delete
 */
void delete_fft_data(struct pvda_fft_data_t* data);

/**
 * @brief Retrieves an NTT table from Prividema-FHE FFT data object
 *
 *
 * @param data The Prividema-FHE FFT data object
 * @param t Modulus for which the table is to be relieved
 *
 * @return A pointer to the vector of NTTRoots if found, NULL otherwise
 */
NTTRoot* get_ntt_table(struct pvda_fft_data_t* data, uint64_t t);

/**
 * @brief Generates and inserts the NTT root table for the specified modulus t
 *
 * @retval 0 on success
 * @retval -1 on failure
 * @retval 1 if t is already present
 *
 */
int generate_ntt_table(const PvdaBackend* backend, uint64_t t);

#endif
