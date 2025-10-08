add_library(COMMON utils.c)
set_target_properties(COMMON PROPERTIES OUTPUT_NAME common)
target_include_directories(COMMON PUBLIC "${PROJECT_SOURCE_DIR}/include/common")
