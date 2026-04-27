# ----------------------------------------
# Build GGSW shared library
# ----------------------------------------

# Source files for GGSW
set(GGSW_SRCS
    ggsw_arithmetic.c
    ggsw_ciphertext.c
    ggsw_params.c
    ggsw_key.c
    glwegadget_arithmetic.c
    glwegadget_ciphertext.c
)

# Create a shared library named GGSW
add_library(PVDA_GGSW SHARED ${GGSW_SRCS})

# Set the output library name explicitly
set_target_properties(PVDA_GGSW PROPERTIES OUTPUT_NAME pvdaggsw)

# Link dependencies: bivGGSW depends on GLWE, libspqlios, and COMMON
target_link_libraries(PVDA_GGSW PRIVATE PVDA_GLWE libspqlios PVDA_COMMON PVDA_BACKEND)

# Include directories
target_include_directories(PVDA_GGSW
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/core/ggsw
)
