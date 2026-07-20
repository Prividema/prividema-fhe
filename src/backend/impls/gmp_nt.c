
#include <stdint.h>

#include "backend.h"
#include "backend_private.h"

uint64_t ref_ring_2nth_root(const PvdaBackend* backend, uint64_t t) { return 0; }

void pvda_fill_gmp_nt(struct pvda_virtual_table* vt) { vt->pvda_ring_2nth_root = ref_ring_2nth_root; }
