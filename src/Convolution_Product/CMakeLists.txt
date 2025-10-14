# Create one shared lib for GLWE
set(CONV_SRCS
    my_vec_znx_fft.c
)

add_library(CONV SHARED ${CONV_SRCS})

set_target_properties(CONV PROPERTIES OUTPUT_NAME conv_product)

target_include_directories(CONV PUBLIC 
        ${PROJECT_SOURCE_DIR}/include/Convolution_Product 
        ${CMAKE_SOURCE_DIR}/build/_deps/sp-src/spqlios/arithmetic)