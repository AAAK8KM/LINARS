message("${CMAKE_CURRENT_SOURCE_DIR}")
# Build local libraries
# ----------------------------------------------

#add_executable( test2 testld.cpp)
#target_link_libraries(test2 SESLOGIC)

include_directories(../include)


if(BUILD_SHARED_LIBS)
    add_library(MATRIX m3diag.cpp)
    target_compile_definitions(MATRIX PUBLIC MATRIX)
    target_compile_features(MATRIX PRIVATE ${TARGET_CXX_STANDARD})
    target_compile_options(MATRIX PRIVATE ${COMPILE_FLAGS})

    add_library(M3DIAG matrix.cpp m3diag.cpp)
    target_compile_definitions(M3DIAG PUBLIC M3DIAG)
    target_compile_features(M3DIAG PRIVATE ${TARGET_CXX_STANDARD})
    target_compile_options(M3DIAG PRIVATE ${COMPILE_FLAGS})
endif()

add_library(MATRIX matrix.cpp m3diag.cpp mdok.cpp mcsr.cpp)
target_compile_definitions(MATRIX PUBLIC MATRIX)
target_compile_features(MATRIX PRIVATE ${TARGET_CXX_STANDARD})
target_compile_options(MATRIX PRIVATE ${COMPILE_FLAGS})

add_library(TOMAS tomas.cpp)
target_compile_definitions(TOMAS PUBLIC TOMAS)
target_link_libraries(TOMAS PRIVATE MATRIX)
target_compile_features(TOMAS PRIVATE ${TARGET_CXX_STANDARD})
target_compile_options(TOMAS PRIVATE ${COMPILE_FLAGS})