
cmake_minimum_required(VERSION 3.14)
project(OMGraphstream)

# file(GLOB OM_GRAPHSTREAM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/gs-netstream/c++/src/*.cpp)

set(OM_GRAPHSTREAM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/gs-netstream/c++/src/netstream-sender.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/gs-netstream/c++/src/netstream-socket.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/gs-netstream/c++/src/netstream-storage.cpp)


add_library(netstream STATIC ${OM_GRAPHSTREAM_SOURCES})
target_link_libraries(netstream INTERFACE wsock32)
target_include_directories(netstream PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/gs-netstream/c++/src/)

add_executable(netstream-main ${CMAKE_CURRENT_SOURCE_DIR}/gs-netstream/c++/src/netstream-main.cpp)
target_link_libraries(netstream-main netstream)
