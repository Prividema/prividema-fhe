cmake_minimum_required(VERSION 3.10)

# Create a library (or change to shared if you want)
#[[
    From what I saw on the top answer by Paolo of this stackoverflow thread https://stackoverflow.com/questions/2649334/difference-between-static-and-shared-libraries :

    - Shared libraries are .so (or in Windows .dll) files. All the code relating to the library is in this file.
    - Static libraries are .a (or in Windows .lib) files. All the code relating to the library is in this file.

    - A shared library is referenced by programs using it at runtime. 
    - A static library is directly linked into the program at compile time.

    - A program using a shared library only makes reference to the code that it uses in the shared library.
    - A program using a static library takes copies of the code that it uses from the static library and makes it part of the program. (Windows also has .lib files which are used to reference .dll files, but they act the same way as the first one).
    
    Advantages in using a shared library :
    - No code duplication => small binary files.
    Using a shared library leads to small binary files, but adds costs during runtime as we have to call the shared library outside of the compiled code.
    Using a static library grants no additionnal runtime loading costs, but generates greater binary files' size.
]]

# can be done similar for the others

# We'll use Shared libraries

# Build core library : The library will be an interface that contains glwe and ggsw.

# Build glwe
add_subdirectory(core/glwe)
add_library(CORE INTERFACE)
target_link_libraries(CORE INTERFACE GLWE)

# Build ggsw
# We define a compiler macro 'GGSW' to avoid errors thrown by codes that need ggsw.
if(BUILD_GGSW)
    target_compile_definitions(CORE INTERFACE GGSW=1)
    add_subdirectory(core/ggsw)
    target_link_libraries(CORE INTERFACE GGSW)
else()
    target_compile_definitions(CORE INTERFACE GGSW=0)
endif(BUILD_GGSW)