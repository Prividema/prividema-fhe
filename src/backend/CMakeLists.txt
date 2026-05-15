# ----------------------------------------
# Build BACKEND library
# ----------------------------------------

# Create BACKEND shared library
add_library(PVDA_BACKEND SHARED
backend.c
spqlios_alias.c
rng.c)

# Set the output library name explicitly
set_target_properties(PVDA_BACKEND PROPERTIES OUTPUT_NAME pvdabackend)

target_link_libraries(PVDA_BACKEND PRIVATE libspqlios PVDA_COMMON OpenSSL::Crypto)

# Include directories
target_include_directories(PVDA_BACKEND
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/backend   # Public headers for consumers
    PRIVATE
        ${PROJECT_SOURCE_DIR}/build/_deps/spqlios-src/spqlios/arithmetic/  # Internal headers only
)
