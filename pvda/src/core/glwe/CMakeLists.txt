# ----------------------------------------
# Build GLWE shared library
# ----------------------------------------

# Source files for GLWE
set(GLWE_SRCS
    glwe_arithmetic.c
    glwe_params.c
    glwe_ciphertext.c
    glwe_key.c
    ${CMAKE_CURRENT_SOURCE_DIR}/../polynomial/univariate_polynomial.c
    ${CMAKE_CURRENT_SOURCE_DIR}/../polynomial/bivariate_polynomial.c
    ${CMAKE_CURRENT_SOURCE_DIR}/../polynomial/random_polynomial.c
)

# Create a shared library named GLWE
add_library(PVDA_GLWE ${GLWE_SRCS})

# Set the output library name explicitly
set_target_properties(PVDA_GLWE PROPERTIES OUTPUT_NAME pvdaglwe)

# Link dependencies
target_link_libraries(PVDA_GLWE PUBLIC PVDA_TYPES PRIVATE PVDA_COMMON PVDA_ADAPTER)

# Include directories
target_include_directories(PVDA_GLWE
    PUBLIC
        ${PVDA_INCLUDE_DIR}/core/glwe/        # Public headers for consumers
        ${PVDA_INCLUDE_DIR}/core/polynomial
)
