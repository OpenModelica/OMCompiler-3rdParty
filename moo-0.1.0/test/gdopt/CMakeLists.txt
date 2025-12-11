cmake_minimum_required(VERSION 3.14)

add_executable(test_gdopt
    entry.cpp
    generated.c
    sanity_check.c
    full_functionality_check.c
)

target_include_directories(test_gdopt PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_gdopt PRIVATE moo)

add_test(NAME test_gdopt COMMAND test_gdopt)
set_tests_properties(test_gdopt PROPERTIES WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
