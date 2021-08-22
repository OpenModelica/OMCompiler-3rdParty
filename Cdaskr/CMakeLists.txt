
cmake_minimum_required(VERSION 3.14)
project(Cdaskr)


set(OM_CDASKR_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/solver/daux.c
                      ${CMAKE_CURRENT_SOURCE_DIR}/solver/ddaskr.c
                      ${CMAKE_CURRENT_SOURCE_DIR}/solver/dlinpk.c)

add_library(cdaskr STATIC)
target_sources(cdaskr PRIVATE ${OM_CDASKR_SOURCES})

install(TARGETS cdaskr)
