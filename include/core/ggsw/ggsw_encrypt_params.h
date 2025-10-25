#ifndef GGSW_ENCRYPT_PARAMS_H
#define GGSW_ENCRYPT_PARAMS_H

#include <stdint.h>

typedef struct ggsw_encrypt_params {
  	int64_t K; 		    // The bivHalfGGSW is a vector of K-normalised and reduced bivRLWE
  	
	int64_t K_tilde;    // each bivHalfGGSW is a family of bivRLWE. 
						// Each bivRLWE encrypts u * 2^(K_tilde*i)
						// For now, same for each bivHalfGGSW
  	
	int64_t L;          // Error precision
	
	uint64_t decomp; 	// Number of limbs per coefficient in a bivRLWE ciphertext 
						// used in an external product
} GGSWEncryptParams;

#endif  // GGSW_ENCRYPT_PARAMS_H
