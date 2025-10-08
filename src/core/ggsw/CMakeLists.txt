# Create one shared lib for GGSW.
set(GGSW_SRCS
    ggsw.c
)
add_library(GGSW SHARED ${GGSW_SRCS})

set_target_properties(GGSW PROPERTIES OUTPUT_NAME ggsw)

# As GGSW requires GLWE we have to link libglwe to the libggsw.
target_link_libraries(GGSW PRIVATE GLWE)

target_include_directories(GGSW PUBLIC ${PROJECT_SOURCE_DIR}/include/core/ggsw)

