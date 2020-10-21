
cmake_minimum_required(VERSION 3.14)

file(GLOB ANTLR_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)


add_library(omantlr3 STATIC ${ANTLR_SOURCES})
target_include_directories(omantlr3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# for antlr3config.h. Needs to be removed after fixing where the file is generated. i.e., move to ./include folder
target_include_directories(omantlr3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/)
target_compile_definitions(omantlr3 PRIVATE ANTLR3_NODEBUGGER
                                            #ANTLR3_USE_64BIT #apparently this is already done by a header file in ./include
                                            )
