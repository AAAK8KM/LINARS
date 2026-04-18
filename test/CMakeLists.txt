
include_directories(../include)

message("${CMAKE_CURRENT_SOURCE_DIR}")

# Build tests
# ----------------------------------------------
add_executable( basicIOTest basic_io_test.cpp)
target_link_libraries(basicIOTest MATRIX)
set_target_properties(basicIOTest PROPERTIES ${TARGET_CXX_PARAM})

add_executable( mulTest mul_test.cpp)
target_link_libraries(mulTest MATRIX)
set_target_properties(mulTest PROPERTIES ${TARGET_CXX_PARAM})

add_executable( tomasTest tomas_test.cpp)
target_link_libraries(tomasTest MATRIX TOMAS)
target_compile_options(tomasTest PRIVATE -O3)
set_target_properties(tomasTest PROPERTIES ${TARGET_CXX_PARAM})


find_package(GTest REQUIRED)
include(GoogleTest)

add_executable( mdokTest mdok_test.cpp)
target_link_libraries(mdokTest PUBLIC MATRIX  PRIVATE GTest::gtest_main)
set_target_properties(mdokTest PROPERTIES ${TARGET_CXX_PARAM})
gtest_discover_tests(mdokTest)

add_executable( mcsrTest mcsr_test.cpp)
target_link_libraries(mcsrTest PUBLIC MATRIX PRIVATE GTest::gtest_main)
set_target_properties(mcsrTest PROPERTIES ${TARGET_CXX_PARAM})
gtest_discover_tests(mcsrTest)

add_executable( vectorTest vector_test.cpp)
target_link_libraries(vectorTest PUBLIC MATRIX PRIVATE GTest::gtest_main)
set_target_properties(vectorTest PROPERTIES ${TARGET_CXX_PARAM})
gtest_discover_tests(vectorTest)

add_executable( mxvTest mxv_test.cpp)
target_link_libraries(mxvTest PUBLIC MATRIX PRIVATE GTest::gtest_main)
set_target_properties(mxvTest PROPERTIES ${TARGET_CXX_PARAM})
gtest_discover_tests(mxvTest)

add_executable( qrTestm qrdec_test.cpp)
target_link_libraries(qrTestm PUBLIC MATRIX PRIVATE GTest::gtest_main)
set_target_properties(qrTestm PROPERTIES ${TARGET_CXX_PARAM})
gtest_discover_tests(qrTestm)

add_executable( qrSolvTest qr_solve_test.cpp)
target_link_libraries(qrSolvTest PUBLIC MATRIX PRIVATE GTest::gtest_main)
set_target_properties(qrSolvTest PROPERTIES ${TARGET_CXX_PARAM})
gtest_discover_tests(qrSolvTest)

add_executable( implictSolvTest implictsolvers_test.cpp)
target_link_libraries(implictSolvTest PUBLIC MATRIX PRIVATE GTest::gtest_main)
set_target_properties(implictSolvTest PROPERTIES ${TARGET_CXX_PARAM})
gtest_discover_tests(implictSolvTest)

add_executable( itTestm it_test_m.cpp)
target_link_libraries(itTestm MATRIX)
set_target_properties(itTestm PROPERTIES ${TARGET_CXX_PARAM})

add_executable( csrTestm csr_test_m.cpp)
target_link_libraries(csrTestm MATRIX)
set_target_properties(csrTestm PROPERTIES ${TARGET_CXX_PARAM})

add_executable( t1Testm time_test.cpp)
target_link_libraries(t1Testm MATRIX)
set_target_properties(t1Testm PROPERTIES ${TARGET_CXX_PARAM})

add_executable( implictTimeTestm time_implict_measeure.cpp)
target_link_libraries(implictTimeTestm MATRIX)
set_target_properties(implictTimeTestm PROPERTIES ${TARGET_CXX_PARAM})

add_executable( implictTimeTestme time_implict_measeure2.cpp)
target_link_libraries(implictTimeTestme MATRIX)
set_target_properties(implictTimeTestme PROPERTIES ${TARGET_CXX_PARAM})

add_executable( ctrl test_1.cpp)
target_link_libraries(ctrl MATRIX)
set_target_properties(ctrl PROPERTIES ${TARGET_CXX_PARAM})

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