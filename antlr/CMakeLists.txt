
cmake_minimum_required(VERSION 3.14)
project(OMAntlr3)

set(OMAntlr3_ANTLRJAR ${CMAKE_CURRENT_SOURCE_DIR}/3.2/tool/antlr-3.2.jar CACHE PATH "Path to antlr jar.")
omc_add_subdirectory(3.2/libantlr3c-3.2)
