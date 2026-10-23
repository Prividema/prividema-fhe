#ifndef SPQLIOS_ALIAS_H
#define SPQLIOS_ALIAS_H

#include <stdint.h>

#include "backend.h"
#include "backend_private.h"

/**
 *
 * @file spqlios_alias.h
 *
 *
 * Functions needed to create and delete a spqlios module, as well as set functions to them in Prividema
 *
 */

SPQLIOS_MODULE* spqlios_new_module_info(uint64_t nn);

void spqlios_delete_module_info(PvdaBackend* module);

void pvda_fill_spqlios(struct pvda_virtual_table* vt);
#endif
