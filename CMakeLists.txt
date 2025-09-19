cmake_minimum_required(VERSION 3.10)
project(prividema-fhe C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

#[[
    Build Options. Watchout for the following point :
    When we build the project and specify the values of the different options, they are registered in a cache at build/CMakeCache.txt
    So if you rebuild the project without specifying any value, thinking it will be OFF by default, Cmake would actualy read these values in its cache.
    For example if I execute :
    ```
        cmake ..  -DBUILD_TESTS=ON
    ```
    and then :
    ```
        cmake ..
    ```
    BUILD_TESTS will still be ON.
]]
option(BUILD_TESTS "Build the Tests" OFF)
option(ENABLE_DEBUG "Enable addtionnal prints for debugging purpose" OFF)

#[[ 
    For setting up debug prints in the code we have 2 choices : 
    On on hand lays Compiler Macros : 
        - A macro is set at compilation and if ENABLE_DEBUG option is off it won't compile the debug code. 
        - But you have to recompile the whole project if you want to turn on the macro.
    On the other hand lays Runtime Variables :
        - We can turn it on/off without recompiling.
        - Tiny performance hit as we have to evaluate all the if(ENABLE_DEBUG).

    For what we do, I think it's better to use compiler macros, even if the code may look a bit messier with lines like '#if ENABLE_DEBUG', '#else' and '#endif'.
]]
if(ENABLE_DEBUG)
    add_compile_definitions(DEBUG=1)
else()
    add_compile_definitions(DEBUG=0)
endif(ENABLE_DEBUG)

# Add include directory
include_directories(include)

# Add subdirectories
add_subdirectory(src)

# Add tests subdirectory if the build option BUILD_TESTS is ON
# The option is OFF by default and can be turned ON by executing cmake with -DBUILD_TESTS=ON
if(BUILD_TESTS)
    add_subdirectory(tests)
endif(BUILD_TESTS)