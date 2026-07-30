#ifndef PVDA_BACKEND_H
#define PVDA_BACKEND_H

#include <stdint.h>

/**
 * @brief Opaque internal type for the prividema backend object type
 */
typedef struct pvda_module_t PvdaBackend;

/**
 * @brief Configuration struct for a Prividema-FHE backend
 *
 * Still very much a workd in progress
 */
typedef struct pvda_module_config_t
{
	uint64_t nn;

} PvdaBackendConfig;

//Legacy backend port functions
//
/**
 * @brief Utility function (and stopgap function until
 * we have the config object system working) to create a
 * backend that uses the spqlios backend
 *
 * @param nn Parameter N
 */
PvdaBackend* pvda_new_spqlios_backend(int nn);

/**
 * NOT YET IMPLEMENTED
 * @brief Create a new backend according to the configuration struct
 *
 *
 * @param config A backend configuration struct
 *
 * @return NULL in case of failure, a pointer to the new backend otherwise
 */
PvdaBackend* pvda_new_backend(PvdaBackendConfig config);

/**
 * @brief Delete a Prividema-FHE backend
 *
 * @param backend The backend to delete
 */
void pvda_delete_backend(PvdaBackend* backend);

/**
 * @brief Retrieves parameter N for a given Prividema-FHE backend
 *
 * @param backend The backend
 *
 * @return Parameter N
 */
uint64_t pvda_module_extract_nn(const PvdaBackend* backend);

#endif
