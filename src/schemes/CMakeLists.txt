# ----------------------------------------
# Build SCHEMES library
# ----------------------------------------

# Create SCHEMES library
add_library(SCHEMES scheme.c)

# Set the output library name explicitly
set_target_properties(SCHEMES PROPERTIES OUTPUT_NAME scheme)

# Include directories
target_include_directories(SCHEMES
    INTERFACE
        ${PROJECT_SOURCE_DIR}/include/schemes   # Public headers for consumers
)
