#include "ggsw_ciphertext.h"

HalfGGSWCiphertext* new_halfggsw(Core* core, uint64_t nrows, uint64_t DegreeInY
){
    HalfGGSWCiphertext* output = malloc(sizeof(HalfGGSWCiphertext*));
    uint64_t N = core->R->nn;
    uint64_t k = core->k;

    output->index = index_halfggsw_vmp;
    output->params = new_halfggsw_encrypt_params(nrows, DegreeInY);
    output->values = new_vmp_pmat(core->R,nrows,DegreeInY*(k+1));
    ((k+1)*DegreeInY*nrows,N*sizeof(int64_t));
    
    return output;
}

GGSWCiphertext* new_ggsw(Core* core, uint64_t K, uint64_t K_tilde, uint64_t L, uint64_t decomp
){
    GGSWCiphertext* output = malloc(sizeof(GGSWCiphertext));
    uint64_t k = core->k;

    // There are k+1 HalfGGSW ciphertexts in a GGSW ciphertext
    output->halfggswVec = malloc((k+1)*sizeof(HalfGGSWCiphertext));

    // Initialize each HalfGGSW ciphertext
    for(uint64_t i = 0 ; i < k+1 ; i++){
        
        output->halfggswVec[i] = new_halfggsw(core, decomp, L/K_tilde);

    }
    // Initialize the parameters of the ggsw ciphertext
    output->params = malloc(sizeof(GGSWEncryptParams));
    
}