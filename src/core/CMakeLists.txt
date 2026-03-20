# ----------------------------------------
# Build GLWE and GGSW modules and CORE interface
# ----------------------------------------

# Build GLWE library
add_subdirectory(glwe)

# Create CORE interface library and link GLWE
add_library(CORE INTERFACE)
target_link_libraries(CORE INTERFACE GLWE)

# Build GGSW library if enabled
# A compiler macro 'GGSW' is defined for code that depends on GGSW
if(BUILD_GGSW)
    target_compile_definitions(CORE INTERFACE GGSW=1)
    add_subdirectory(ggsw)
    target_link_libraries(CORE INTERFACE GGSW)
else()
    target_compile_definitions(CORE INTERFACE GGSW=0)
endif()
