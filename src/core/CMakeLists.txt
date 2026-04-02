# ----------------------------------------
# Build GLWE and GGSW modules and CORE interface
# ----------------------------------------

# Build GLWE library
add_subdirectory(glwe)
add_subdirectory(ggsw)

# Create CORE interface library and link GLWE
add_library(CORE INTERFACE)
target_link_libraries(CORE INTERFACE GLWE GGSW)

