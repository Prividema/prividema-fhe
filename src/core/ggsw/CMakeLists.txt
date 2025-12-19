# Create one shared lib for GGSW.
set(GGSW_SRCS
    ggsw.c
    ggsw_ciphertext.c
    ggsw_ct_params.c
)
add_library(GGSW SHARED ${GGSW_SRCS})

set_target_properties(GGSW PROPERTIES OUTPUT_NAME ggsw)

# As GGSW requires GLWE we have to link libglwe to the libggsw.
target_link_libraries(GGSW PRIVATE GLWE libspqlios COMMON)

target_include_directories(GGSW PUBLIC  ${PROJECT_SOURCE_DIR}/include/core/ggsw 
                                        ${PROJECT_SOURCE_DIR}/include/common/
                                PRIVATE ${PROJECT_SOURCE_DIR}/build/_deps/spqlios-src/spqlios/arithmetic/)

