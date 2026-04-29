# ----------------------------------------
# Build COMMON shared library
# ----------------------------------------

# Source files for COMMON
set(COMMON_SRCS
    logger.c
    maths_structures.c
    utils.c
)

# Create a shared library named COMMON
add_library(PVDA_COMMON SHARED ${COMMON_SRCS})

# Set the output library name explicitly
set_target_properties(PVDA_COMMON PROPERTIES OUTPUT_NAME pvdacommon)

# Link dependencies
target_link_libraries(PVDA_COMMON PRIVATE)

# Include directories
target_include_directories(PVDA_COMMON
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/common                 # Public headers for consumers
)
