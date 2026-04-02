# ----------------------------------------
# Build GLWE and GGSW modules and CORE interface
# ----------------------------------------

# Build GLWE library
add_subdirectory(glwe)
add_subdirectory(ggsw)

# Create CORE interface library and link GLWE
add_library(PVDA_CORE INTERFACE)
target_link_libraries(PVDA_CORE INTERFACE PVDA_GLWE PVDA_GGSW PVDA_COMMON)

