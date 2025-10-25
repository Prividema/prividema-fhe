#ifndef HALFGLWE_ENCRYPT_PARAMS_H
#define HALFGLWE_ENCRYPT_PARAMS_H

#include <stdint.h>

/** @brief Computes the index of the 
 *         jjj-ème coefficient of the 
 *         jj-ème poly of the 
 *         j-ème bivariate poly of the 
 *         i-ème bivRLWE  */

typedef struct halfggsw_encrypt_params{
    uint64_t nrows;         // Number of term 
    uint64_t DegreeInY;     // Maximum degree in Y
} HalfGGSWEncryptParams;

HalfGGSWEncryptParams* new_halfggsw_encrypt_params( uint64_t nrows, uint64_t DegreeInY);

#endif  // HALFGLWE_ENCRYPT_PARAMS_H