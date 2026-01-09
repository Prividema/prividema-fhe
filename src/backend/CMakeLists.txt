# ----------------------------------------
# Build BACKEND library
# ----------------------------------------

# Create BACKEND shared library
add_library(BACKEND SHARED backend.c)

# Set the output library name explicitly
set_target_properties(BACKEND PROPERTIES OUTPUT_NAME backend)

# Include directories
target_include_directories(BACKEND
    INTERFACE
        ${PROJECT_SOURCE_DIR}/include/backend   # Public headers for consumers
)
