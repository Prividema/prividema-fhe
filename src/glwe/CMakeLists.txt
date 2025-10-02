# Create one shared lib for GLWE
file(GLOB GLWE_SRCS "*.c")
add_library(GLWE SHARED ${GLWE_SRCS})
target_include_directories(GLWE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include/glwe)