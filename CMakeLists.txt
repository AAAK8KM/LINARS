cmake_minimum_required(VERSION 3.7)
project( LINARS  LANGUAGES CXX)
set (CMAKE_GENERATOR "Unix Makefiles" CACHE INTERNAL "" FORCE)

# Set standarts
# ----------------------------------------------
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic -O2)

# got to the code stuff
include(./src/CMakeLists.txt)
include(./test/CMakeLists.txt)

add_custom_target(mtest ALL
    COMMAND ${CMAKE_CTEST_COMMAND}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running CTest suite"
    VERBATIM 
)
add_dependencies(mtest ${TDEPS})