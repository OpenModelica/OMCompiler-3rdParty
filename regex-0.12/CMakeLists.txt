cmake_minimum_required(VERSION 3.14)

add_library(omcregex STATIC regex.c)
target_include_directories(omcregex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
