message("${CMAKE_CURRENT_SOURCE_DIR}")
# Build local libraries
# ----------------------------------------------

#add_executable( test2 ${CMAKE_CURRENT_SOURCE_DIR}/testld.cpp)
#target_link_libraries(test2 SESLOGIC)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(MATRIX ${CMAKE_CURRENT_SOURCE_DIR}/matrix.cpp)
target_compile_definitions(MATRIX PUBLIC MATRIX)
target_compile_features(MATRIX PRIVATE ${TARGET_CXX_STANDARD})
target_compile_options(MATRIX PRIVATE ${COMPILE_FLAGS})

add_library(TOMAS ${CMAKE_CURRENT_SOURCE_DIR}/tomas.cpp)
target_compile_definitions(TOMAS PUBLIC TOMAS)
target_link_libraries(TOMAS PRIVATE MATRIX)
target_compile_features(TOMAS PRIVATE ${TARGET_CXX_STANDARD})
target_compile_options(TOMAS PRIVATE ${COMPILE_FLAGS})