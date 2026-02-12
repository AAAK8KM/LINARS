set(TMP ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)


message("${CMAKE_CURRENT_SOURCE_DIR}")
# Build local libraries
# ----------------------------------------------

add_executable( basicIOTest ${CMAKE_CURRENT_SOURCE_DIR}/basic_io_test.cpp)
target_link_libraries(basicIOTest MATRIX)

add_executable( mulTest ${CMAKE_CURRENT_SOURCE_DIR}/mul_test.cpp)
target_link_libraries(mulTest MATRIX)


add_executable( tomasTest ${CMAKE_CURRENT_SOURCE_DIR}/tomas_test.cpp)
target_link_libraries(tomasTest MATRIX TOMAS)

set(TDEPS basicIOTest mulTest tomasTest)

set(TEST_CONFIGS
    "basicIOTest d ${CMAKE_CURRENT_SOURCE_DIR}/test/io1"
    "basicIOTest d ${CMAKE_CURRENT_SOURCE_DIR}/test/io2"
    "basicIOTest d ${CMAKE_CURRENT_SOURCE_DIR}/test/io3"
    "basicIOTest d ${CMAKE_CURRENT_SOURCE_DIR}/test/io4"
    "basicIOTest f ${CMAKE_CURRENT_SOURCE_DIR}/test/io1"
    "basicIOTest f ${CMAKE_CURRENT_SOURCE_DIR}/test/io2"
    "basicIOTest f ${CMAKE_CURRENT_SOURCE_DIR}/test/io3"
    "basicIOTest f ${CMAKE_CURRENT_SOURCE_DIR}/test/io4"
    "mulTest d ${CMAKE_CURRENT_SOURCE_DIR}/test/mul1"
    "mulTest d ${CMAKE_CURRENT_SOURCE_DIR}/test/mul2"
    "mulTest d ${CMAKE_CURRENT_SOURCE_DIR}/test/mul3"
    "mulTest f ${CMAKE_CURRENT_SOURCE_DIR}/test/mul1"
    "mulTest f ${CMAKE_CURRENT_SOURCE_DIR}/test/mul2"
    "mulTest f ${CMAKE_CURRENT_SOURCE_DIR}/test/mul3"
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

set(CMAKE_CURRENT_SOURCE_DIR  ${TMP})
