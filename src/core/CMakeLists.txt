
# Build glwe
add_subdirectory(glwe)
add_library(CORE INTERFACE)
target_link_libraries(CORE INTERFACE GLWE)

# Build ggsw
# We define a compiler macro 'GGSW' to avoid errors thrown by codes that need ggsw.
if(BUILD_GGSW)
    target_compile_definitions(CORE INTERFACE GGSW=1)
    add_subdirectory(ggsw)
    target_link_libraries(CORE INTERFACE GGSW)
else()
    target_compile_definitions(CORE INTERFACE GGSW=0)
endif(BUILD_GGSW)