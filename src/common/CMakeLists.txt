# ----------------------------------------
# Build COMMON shared library
# ----------------------------------------

# Source files for COMMON
set(COMMON_SRCS 
    utils.c 
    distributions.c
    structure_alias.c
    spqlios_alias.c
)

# Create a shared library named COMMON
add_library(COMMON SHARED ${COMMON_SRCS})

# Set the output library name explicitly
set_target_properties(COMMON PROPERTIES OUTPUT_NAME common)

# Link dependencies
target_link_libraries(COMMON PRIVATE libspqlios)

# Include directories
target_include_directories(COMMON
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/common                 # Public headers for consumers
    PRIVATE
        ${PROJECT_SOURCE_DIR}/build/_deps/spqlios-src/spqlios/arithmetic/  # Internal headers only
)
