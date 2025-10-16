cmake_minimum_required(VERSION 3.10)

# Build backend library : 
add_subdirectory(backend)

# Build common library :
add_subdirectory(common)

# Build core library : The library will be an interface that contains glwe and ggsw.
add_subdirectory(core)

# Build schemes library :
add_subdirectory(schemes)
