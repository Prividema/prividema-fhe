set(COMMON_SRCS 
    utils.c 
    distributions.c
    spqlios_alias.c)
add_library(COMMON ${COMMON_SRCS})
set_target_properties(COMMON PROPERTIES OUTPUT_NAME common)
target_include_directories(COMMON PUBLIC "${PROJECT_SOURCE_DIR}/include/common")