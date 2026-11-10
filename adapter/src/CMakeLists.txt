# ----------------------------------------
# Build ADAPTER library
# ----------------------------------------

add_library(PVDA_ADAPTER
backend.c
impls/spqlios_alias.c
impls/cpu_nt.c
rng.c
impls/rng_cpu.c
pvda_ffts.c)

# Set the output library name explicitly
set_target_properties(PVDA_ADAPTER PROPERTIES OUTPUT_NAME pvdaadapter)

# The adapter API is expressed with prividema types
target_link_libraries(PVDA_ADAPTER PUBLIC PVDA_TYPES PRIVATE PVDA_COMMON)

if(BUILD_SHARED_LIBS)
  target_link_libraries(PVDA_ADAPTER PRIVATE libspqlios)
else()
  target_link_libraries(PVDA_ADAPTER PRIVATE libspqlios-static)
endif()

if(OPENSSL_BACKEND)

  find_package(OpenSSL REQUIRED)

  target_sources(PVDA_ADAPTER PRIVATE impls/rng_openssl.c)

  target_link_libraries(PVDA_ADAPTER PRIVATE OpenSSL::Crypto)

  target_compile_definitions(PVDA_ADAPTER PRIVATE OPENSSL_BACKEND)

endif()

# Include directories
target_include_directories(PVDA_ADAPTER
    PUBLIC
        ${ADAPTER_INCLUDE_DIR}   # Public headers for consumers
    PRIVATE
        ${PROJECT_SOURCE_DIR}/backends/spqlios-arithmetic/spqlios/arithmetic/  # Internal headers only
)
