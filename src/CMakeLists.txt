cmake_minimum_required(VERSION 3.10)

# Collect all .c files in src
file(GLOB SRC_FILES *.c)

# Create a library (or change to shared if you want)
add_library(core ${SRC_FILES})