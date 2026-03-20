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
    bivariate_polynomial.c
)

# Create a shared library named GLWE
add_library(GLWE SHARED ${GLWE_SRCS})

# Set the output library name explicitly
set_target_properties(GLWE PROPERTIES OUTPUT_NAME glwe)

# Link dependencies
target_link_libraries(GLWE PRIVATE libspqlios COMMON)

# Include directories
target_include_directories(GLWE
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/core/glwe/   # Public headers for consumers
        ${PROJECT_SOURCE_DIR}/include/common/
    PRIVATE
        ${PROJECT_SOURCE_DIR}/build/_deps/spqlios-src/spqlios/arithmetic/  # Internal headers only
)
