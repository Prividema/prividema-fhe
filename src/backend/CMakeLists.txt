# SPQLios
add_library(BACKEND SHARED backend.c)
set_target_properties(BACKEND PROPERTIES OUTPUT_NAME backend)
target_include_directories(BACKEND INTERFACE "${PROJECT_SOURCE_DIR}/include/backend")
