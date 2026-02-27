message("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(../include)


# Build local libraries
# ----------------------------------------------
if(BUILD_SHARED_LIBS)
    add_library(DMATRIX matrix.cpp)
    target_compile_definitions(DMATRIX PUBLIC DMATRIX)
    target_compile_features(DMATRIX PRIVATE ${TARGET_CXX_STANDARD})
    target_compile_options(DMATRIX PRIVATE ${COMPILE_FLAGS})

    add_library(M3DIAG m3diag.cpp)
    target_compile_definitions(M3DIAG PUBLIC M3DIAG)
    target_compile_features(M3DIAG PRIVATE ${TARGET_CXX_STANDARD})
    target_compile_options(M3DIAG PRIVATE ${COMPILE_FLAGS})

    add_library(MDOK mdok.cpp)
    target_compile_definitions(MDOK PUBLIC MDOK)
    target_compile_features(MDOK PRIVATE ${TARGET_CXX_STANDARD})
    target_compile_options(MDOK PRIVATE ${COMPILE_FLAGS})

    add_library(MCSR m3diag.cpp)
    target_compile_definitions(MCSR PUBLIC MCSR)
    target_compile_features(MCSR PRIVATE ${TARGET_CXX_STANDARD})
    target_compile_options(MCSR PRIVATE ${COMPILE_FLAGS})
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