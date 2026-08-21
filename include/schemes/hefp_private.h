#ifndef HEFP_PRIVATE
#define HEFP_PRIVATE

#include "backend.h"
void hefp_encode_slow_internal(double* out, uint64_t n, double _Complex* in);
void hefp_decode_slow_internal(double _Complex* out, uint64_t n, double* in);
int hefp_encode_internal(const PvdaBackend* backend, double* out, uint64_t slots, double _Complex* in, int inplace);
int hefp_decode_internal(const PvdaBackend* backend, double _Complex* out, uint64_t slots, double* in);
int hefp_decode_fft(const PvdaBackend* backend, double _Complex* inout, uint64_t slots);
void hefp_decode_slow_internal(double _Complex* out, uint64_t n, double* in);
#endif
