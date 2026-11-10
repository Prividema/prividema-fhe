# ----------------------------------------
# PVDA: the prividema library
# ----------------------------------------

# Types owned by prividema (polynomials and opaque BE_* backend handles).
# Header-only, shared with the adapter.
add_library(PVDA_TYPES INTERFACE)
target_include_directories(PVDA_TYPES INTERFACE ${PVDA_INCLUDE_DIR}/types)

add_subdirectory(src)

if(BUILD_TESTS)
  add_subdirectory(tests)
endif()
