
include_directories(../include)

message("${CMAKE_CURRENT_SOURCE_DIR}")

# Build tests
# ----------------------------------------------
add_executable( basicIOTest basic_io_test.cpp)
target_link_libraries(basicIOTest MATRIX)
target_compile_features(basicIOTest PRIVATE ${TARGET_CXX_STANDARD})

add_executable( mulTest mul_test.cpp)
target_link_libraries(mulTest MATRIX)
target_compile_features(mulTest PRIVATE ${TARGET_CXX_STANDARD})

add_executable( tomasTest tomas_test.cpp)
target_link_libraries(tomasTest MATRIX TOMAS)
target_compile_options(tomasTest PRIVATE -O3)
target_compile_features(tomasTest PRIVATE ${TARGET_CXX_STANDARD})


find_package(GTest REQUIRED)
include(GoogleTest)

add_executable( mdokTest mdok_test.cpp)
target_link_libraries(mdokTest PUBLIC MATRIX  PRIVATE GTest::gtest_main)
target_compile_features(mdokTest PRIVATE ${TARGET_CXX_STANDARD})
gtest_discover_tests(mdokTest)

add_executable( mcsrTest mcsr_test.cpp)
target_link_libraries(mcsrTest PUBLIC MATRIX PRIVATE GTest::gtest_main)
target_compile_features(mcsrTest PRIVATE ${TARGET_CXX_STANDARD})
gtest_discover_tests(mcsrTest)

add_executable( vectorTest vector_test.cpp)
target_link_libraries(vectorTest PUBLIC MATRIX PRIVATE GTest::gtest_main)
target_compile_features(vectorTest PRIVATE ${TARGET_CXX_STANDARD})
gtest_discover_tests(vectorTest)

add_executable( mxvTest mxv_test.cpp)
target_link_libraries(mxvTest PUBLIC MATRIX PRIVATE GTest::gtest_main)
target_compile_features(mxvTest PRIVATE ${TARGET_CXX_STANDARD})
gtest_discover_tests(mxvTest)

add_executable( qrTestm qrdec_test.cpp)
target_link_libraries(qrTestm PUBLIC MATRIX PRIVATE GTest::gtest_main)
target_compile_features(qrTestm PRIVATE ${TARGET_CXX_STANDARD})
gtest_discover_tests(qrTestm)

add_executable( qrSolvTest qr_solve_test.cpp)
target_link_libraries(qrSolvTest PUBLIC MATRIX PRIVATE GTest::gtest_main)
target_compile_features(qrSolvTest PRIVATE ${TARGET_CXX_STANDARD})
gtest_discover_tests(qrSolvTest)

add_executable( implictSolvTest implictsolvers_test.cpp)
target_link_libraries(implictSolvTest PUBLIC MATRIX PRIVATE GTest::gtest_main)
target_compile_features(implictSolvTest PRIVATE ${TARGET_CXX_STANDARD})
gtest_discover_tests(implictSolvTest)

add_executable( itTestm it_test_m.cpp)
target_link_libraries(itTestm MATRIX)
target_compile_features(itTestm PRIVATE ${TARGET_CXX_STANDARD})

add_executable( csrTestm csr_test_m.cpp)
target_link_libraries(csrTestm MATRIX)
target_compile_features(csrTestm PRIVATE ${TARGET_CXX_STANDARD})

add_executable( t1Testm time_test.cpp)
target_link_libraries(t1Testm MATRIX)
target_compile_features(t1Testm PRIVATE ${TARGET_CXX_STANDARD})

add_executable( implictTimeTestm time_implict_measeure.cpp)
target_link_libraries(implictTimeTestm MATRIX)
target_compile_features(implictTimeTestm PRIVATE ${TARGET_CXX_STANDARD})

set(TDEPS basicIOTest mulTest tomasTest mcsrTest mdokTest mxvTest vectorTest qrSolvTest qrTestm implictSolvTest)

set(TEST_CONFIGS
    "basicIOTest d test/io1"
    "basicIOTest d test/io2"
    "basicIOTest d test/io3"
    "basicIOTest d test/io4"
    "basicIOTest f test/io1"
    "basicIOTest f test/io2"
    "basicIOTest f test/io3"
    "basicIOTest f test/io4"
    "mulTest d test/mul1"
    "mulTest d test/mul2"
    "mulTest d test/mul3"
    "mulTest f test/mul1"
    "mulTest f test/mul2"
    "mulTest f test/mul3"
    "tomasTest d 10 10"
    "tomasTest d 100 100"
    "tomasTest d 1000 100"
    "tomasTest d 10000 2"
    "tomasTest ld 10000 2"
    "tomasTest ld 50000 1"
)

foreach(config IN LISTS TEST_CONFIGS)
    # Split config: first element = test name, rest = args
    
    list(FIND TEST_CONFIGS ${config} idx)
    string(REPLACE " " ";" config ${config})
    list(GET config 0 test_name)
    list(REMOVE_AT config 0)
    
    message("test ${test_name}-${idx} added")
    add_test(
        NAME "${test_name}-${idx}"
        COMMAND ${test_name} ${config}
    )
endforeach()

add_custom_target(mtest
    COMMAND ${CMAKE_CTEST_COMMAND} -j ${CMAKE_BUILD_PARALLEL_LEVEL}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
    COMMENT "Running CTest suite"
    VERBATIM 
)

add_dependencies(mtest ${TDEPS})