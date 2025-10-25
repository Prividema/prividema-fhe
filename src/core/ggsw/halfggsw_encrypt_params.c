#include "halfggsw_encrypt_params.h"

HalfGGSWEncryptParams* new_halfggsw_encrypt_params( uint64_t nrows, uint64_t DegreeInY
){
    HalfGGSWEncryptParams* output = malloc(sizeof(HalfGGSWEncryptParams));
    
    output->nrows = nrows;
    output->DegreeInY = DegreeInY;

    return output;
}