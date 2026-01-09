# ----------------------------------------
# Build GGSW shared library
# ----------------------------------------

# Source files for GGSW
set(GGSW_SRCS
    ggsw.c
    ggsw_ciphertext.c
    ggsw_ct_params.c
)

# Create a shared library named GGSW
add_library(GGSW SHARED ${GGSW_SRCS})

# Set the output library name explicitly
set_target_properties(GGSW PROPERTIES OUTPUT_NAME ggsw)

# Link dependencies: GGSW depends on GLWE, libspqlios, and COMMON
target_link_libraries(GGSW PRIVATE GLWE libspqlios COMMON)

# Include directories
target_include_directories(GGSW
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/core/ggsw    # Public headers for consumers
        ${PROJECT_SOURCE_DIR}/include/common/
    PRIVATE
        ${PROJECT_SOURCE_DIR}/build/_deps/spqlios-src/spqlios/arithmetic/  # Internal headers only
)
