# ----------------------------------------
# Build BACKEND library
# ----------------------------------------

# Create BACKEND shared library
add_library(PVDA_BACKEND
backend.c
impls/spqlios_alias.c
impls/cpu_nt.c
rng.c
impls/rng_cpu.c
pvda_ffts.c)

# Set the output library name explicitly
set_target_properties(PVDA_BACKEND PROPERTIES OUTPUT_NAME pvdabackend)

target_link_libraries(PVDA_BACKEND PRIVATE PVDA_COMMON)

if(BUILD_SHARED_LIBS)
  target_link_libraries(PVDA_BACKEND PRIVATE libspqlios)
else()
  target_link_libraries(PVDA_BACKEND PRIVATE libspqlios-static)
endif()



if(OPENSSL_BACKEND)

  find_package(OpenSSL REQUIRED)

  target_sources(PVDA_BACKEND PRIVATE impls/rng_openssl.c)

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
