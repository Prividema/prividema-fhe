# ----------------------------------------
# Build GLWE shared library
# ----------------------------------------

# Source files for GLWE
set(GLWE_SRCS
    glwe.c
    glwe_params.c
    glwe_ciphertext.c
    glwe_key.c
    glwe_transform_key.c
    univariate_polynomial.c
    bivariate_polynomial.c
)

# Create a shared library named GLWE
add_library(PVDA_GLWE SHARED ${GLWE_SRCS})

# Set the output library name explicitly
set_target_properties(PVDA_GLWE PROPERTIES OUTPUT_NAME pvdaglwe)

# Link dependencies
target_link_libraries(PVDA_GLWE PRIVATE libspqlios PVDA_COMMON)

# Include directories
target_include_directories(PVDA_GLWE
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/core/glwe/   # Public headers for consumers
)
