add_library(SCHEMES scheme.c)
set_target_properties(SCHEMES PROPERTIES OUTPUT_NAME scheme)
target_include_directories(SCHEMES INTERFACE "${PROJECT_SOURCE_DIR}/include/schemes")
