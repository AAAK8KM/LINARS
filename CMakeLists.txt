cmake_minimum_required(VERSION 3.45)
project( LINARS  LANGUAGES CXX)

# Set standarts
# ----------------------------------------------
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(TARGET_CXX_STANDARD cxx_std_23)
if(MSVC)
    set(COMPILE_FLAGS /W4 /WX /permissive-)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang|Intel")
    set(COMPILE_FLAGS -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Werror)
endif()

# got to the code stuff
add_subdirectory(src)
option(TESTING "enable testing" ON)
if (TESTING)
    enable_testing()
    add_subdirectory(test)
endif()