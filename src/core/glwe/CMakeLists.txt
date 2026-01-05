# Create one shared lib for GLWE
set(GLWE_SRCS
    glwe.c
    glwe_ct_params.c
    glwe_ciphertext.c
)

add_library(GLWE SHARED ${GLWE_SRCS})

set_target_properties(GLWE PROPERTIES OUTPUT_NAME glwe)

target_link_libraries(GLWE PRIVATE libspqlios COMMON)
target_include_directories(GLWE PUBLIC  ${PROJECT_SOURCE_DIR}/include/core/glwe/
                                        ${PROJECT_SOURCE_DIR}/include/common/
                                PRIVATE ${PROJECT_SOURCE_DIR}/build/_deps/spqlios-src/spqlios/arithmetic/)