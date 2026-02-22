
include_directories(../include)

#set(CMAKE_CURRENT_SOURCE_DIR test)


message("${CMAKE_CURRENT_SOURCE_DIR}")
# Build local libraries
# ----------------------------------------------

add_executable( basicIOTest basic_io_test.cpp)
target_link_libraries(basicIOTest MATRIX)
target_compile_features(basicIOTest PRIVATE ${TARGET_CXX_STANDARD})

add_executable( mulTest mul_test.cpp)
target_link_libraries(mulTest MATRIX)
target_compile_features(mulTest PRIVATE ${TARGET_CXX_STANDARD})

add_executable( tomasTest tomas_test.cpp)
target_link_libraries(tomasTest MATRIX TOMAS)
target_compile_features(tomasTest PRIVATE ${TARGET_CXX_STANDARD})

set(TDEPS basicIOTest mulTest tomasTest)

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

enable_testing()

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

add_dependencies(mtest ${TDEPS})