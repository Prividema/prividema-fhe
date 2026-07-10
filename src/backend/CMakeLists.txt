# ----------------------------------------
# Build BACKEND library
# ----------------------------------------

# Create BACKEND shared library
add_library(PVDA_BACKEND
backend.c
spqlios_alias.c
rng.c
rng_cpu.c
pvda_ffts.c)

# Set the output library name explicitly
set_target_properties(PVDA_BACKEND PROPERTIES OUTPUT_NAME pvdabackend)

target_link_libraries(PVDA_BACKEND PRIVATE libspqlios PVDA_COMMON)

if(OPENSSL_BACKEND)

  find_package(OpenSSL REQUIRED)

  target_sources(PVDA_BACKEND PRIVATE rng_openssl.c)

  target_link_libraries(PVDA_BACKEND PRIVATE  OpenSSL::Crypto)

  target_compile_definitions(PVDA_BACKEND PRIVATE OPENSSL_BACKEND)

endif()

# Include directories
target_include_directories(PVDA_BACKEND
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/backend   # Public headers for consumers
    PRIVATE
        ${PROJECT_SOURCE_DIR}/spqlios-arithmetic/spqlios/arithmetic/  # Internal headers only
)
