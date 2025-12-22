set(COMMON_SRCS 
    utils.c 
    distributions.c
    structure_alias.c
    spqlios_alias.c
)
add_library(COMMON SHARED ${COMMON_SRCS})
set_target_properties(COMMON PROPERTIES OUTPUT_NAME common)
target_link_libraries(COMMON PRIVATE libspqlios)
target_include_directories(COMMON PUBLIC ${PROJECT_SOURCE_DIR}/include/common
                                         ${PROJECT_SOURCE_DIR}/build/_deps/spqlios-src/spqlios/arithmetic/)