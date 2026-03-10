#include "ggsw_key.h"

#include "glwe_key.h"
#include "logger.h"
#include "rng.h"
#include "utils.h"
#include "spqlios_alias.h"

//! bivGGSW PART (begin)

PolyUniv** new_ggsw_secret_key_values(uint64_t N, uint64_t k)
{
  return new_glwe_secret_key_values(N, k);
}

void delete_ggsw_secret_key_values(PolyUniv** values, uint64_t k)
{
  delete_glwe_secret_key_values(values, k);
}

GGSWSecretKey* new_ggsw_secret_key(uint64_t N, uint64_t k)
{
  return new_glwe_secret_key(N, k);
}

int uniform_ggsw_secret_key(const MODULE* module, GGSWSecretKey* sk, uint64_t nb_bits)
{
  return uniform_glwe_secret_key(module, sk, nb_bits);
}

void delete_ggsw_secret_key(GGSWSecretKey* sk)
{
  delete_glwe_secret_key(sk);
}


PolyUnivDFT** new_ggsw_secret_key_values_dft(uint64_t N, uint64_t k)
{
  return new_glwe_secret_key_values_dft(N,k);
}

void delete_ggsw_secret_key_values_dft(PolyUnivDFT** values, uint64_t k)
{
  delete_glwe_secret_key_values_dft(values, k);
}

GGSWSecretKeyDFT* new_ggsw_secret_key_dft(uint64_t N, uint64_t k)
{
  return new_glwe_secret_key_dft(N, k);
}

int uniform_ggsw_secret_key_dft(const MODULE* module, GGSWSecretKeyDFT* sk_dft, uint64_t nb_bits)
{
  return uniform_glwe_secret_key_dft(module, sk_dft, nb_bits);
}

void delete_ggsw_secret_key_dft(GGSWSecretKeyDFT* sk_dft)
{
  delete_glwe_secret_key_dft(sk_dft);
}
