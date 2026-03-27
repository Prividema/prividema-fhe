# ----------------------------------------
# Build GGSW shared library
# ----------------------------------------

# Source files for GGSW
set(GGSW_SRCS
    ggsw.c
    ggsw_ciphertext.c
    ggsw_params.c
    ggsw_key.c
    glwegadget.c
    glwegadget_ciphertext.c
)

# Create a shared library named GGSW
add_library(GGSW SHARED ${GGSW_SRCS})

# Set the output library name explicitly
set_target_properties(GGSW PROPERTIES OUTPUT_NAME ggsw)

# Link dependencies: bivGGSW depends on GLWE, libspqlios, and COMMON
target_link_libraries(GGSW PRIVATE GLWE libspqlios COMMON)

# Include directories
target_include_directories(GGSW
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/core/ggsw
)
