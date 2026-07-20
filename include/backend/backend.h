#ifndef PVDA_BACKEND_H
#define PVDA_BACKEND_H

#include <stdint.h>
#include <sys/types.h>

#include "maths_structures.h"

typedef struct pvda_module_t PvdaBackend;

typedef struct pvda_module_config_t
{
	uint64_t nn;

} PvdaBackendConfig;

//Legacy backend port functions
PvdaBackend* pvda_new_spqlios_backend(int nn);

void pvda_delete_backend(PvdaBackend* backend);

uint64_t pvda_module_extract_nn(const PvdaBackend* backend);

#endif
