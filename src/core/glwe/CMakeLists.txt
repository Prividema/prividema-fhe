# Create one shared lib for GLWE
set(GLWE_SRCS
    glwe.c
    glwe_ct_params.c
)

add_library(GLWE SHARED ${GLWE_SRCS})

set_target_properties(GLWE PROPERTIES OUTPUT_NAME glwe)

target_include_directories(GLWE PUBLIC ${PROJECT_SOURCE_DIR}/include/core/glwe)