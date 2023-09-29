cmake_minimum_required(VERSION 3.14)
project(ryu)

add_library(ryu STATIC)
target_sources(ryu PRIVATE ryu/d2fixed.c  ryu/d2s.c ryu/f2s.c ryu/generic_128.c ryu/s2d.c ryu/s2f.c ryu/om_format.c)

target_include_directories(ryu PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/ryu)