cmake_minimum_required(VERSION 3.14)
project(opcua)

add_library(opcua STATIC)
target_sources(opcua PRIVATE open62541.c)

target_include_directories(opcua PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
