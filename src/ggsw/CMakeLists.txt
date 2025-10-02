# Create one shared lib for GGSW.
file(GLOB GGSW_SRCS "*.c")
add_library(GGSW SHARED ${GGSW_SRCS})
target_include_directories(GGSW PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include/ggsw
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include/glwe
)

# As GGSW requires GLWE we have to link libGLWE to the libGGSW.
target_link_libraries(GGSW PUBLIC GLWE)