# ----------------------------------------
# Build SCHEMES library
# ----------------------------------------

# Create SCHEMES library
add_library(PVDA_SCHEMES
scheme.c
tfhe.c
hefp.c
hefp_internals.c)

# Set the output library name explicitly
set_target_properties(PVDA_SCHEMES PROPERTIES OUTPUT_NAME pvda_schemes)

target_link_libraries(PVDA_SCHEMES PRIVATE libspqlios PVDA_COMMON PVDA_BACKEND PVDA_CORE)
# Include directories
target_include_directories(PVDA_SCHEMES
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/schemes   # Public headers for consumers
)
