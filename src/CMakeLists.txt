set(TMP ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

message("${CMAKE_CURRENT_SOURCE_DIR}")
# Build local libraries
# ----------------------------------------------

#add_executable( test2 ${CMAKE_CURRENT_SOURCE_DIR}/testld.cpp)
#target_link_libraries(test2 SESLOGIC)

add_library(MATRIX SHARED ${CMAKE_CURRENT_SOURCE_DIR}/matrix.cpp)
target_compile_definitions(MATRIX PUBLIC MATRIX)


add_library(TOMAS SHARED ${CMAKE_CURRENT_SOURCE_DIR}/tomas.cpp)
target_compile_definitions(TOMAS PUBLIC TOMAS)
target_link_libraries(TOMAS PUBLIC MATRIX)



set(CMAKE_CURRENT_SOURCE_DIR  ${TMP})
