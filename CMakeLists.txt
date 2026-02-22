cmake_minimum_required(VERSION 3.7)
project( LINARS  LANGUAGES CXX)

# Set standarts
# ----------------------------------------------
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(TARGET_CXX_STANDARD cxx_std_20)
#Additional compile flags only for this project
if(MSVC)
        set(COMPILE_FLAGS /W4 /WX /permissive-)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang|Intel")
    set(COMPILE_FLAGS -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Werror)
endif()

add_custom_target(mtest ALL
    COMMAND ${CMAKE_CTEST_COMMAND}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
    COMMENT "Running CTest suite"
    VERBATIM 
)

# got to the code stuff

add_subdirectory(src)
add_subdirectory(test)