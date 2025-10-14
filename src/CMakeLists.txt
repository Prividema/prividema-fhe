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
# Build backend library : 
add_subdirectory(backend)

# Build common library :
add_subdirectory(common)

# Build Conv_Prod library
# Build core library : The library will be an interface that contains glwe and ggsw.
add_subdirectory(core)

# Build schemes library :
add_subdirectory(schemes)

add_subdirectory(Convolution_Product)

include(FetchContent)

FetchContent_Declare(
  sp
  GIT_REPOSITORY https://github.com/tfhe/spqlios-arithmetic.git
  # Choisir une branche ou un tag existant
  GIT_TAG main
)

FetchContent_MakeAvailable(sp)