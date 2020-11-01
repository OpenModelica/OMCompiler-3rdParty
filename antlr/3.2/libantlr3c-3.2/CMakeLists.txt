
cmake_minimum_required(VERSION 3.14)

file(GLOB ANTLR_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)

option(ANTLR3_NODEBUGGER "Turns off default flags that include the antlr debugger in the runtime. Specify to remove debugger and the socket dependancies" ON)
option(ANTLR3_USE_64BIT "Turns on flags that produce 64 bit object code if any are required" ON)

include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckIncludeFile)
include(CheckTypeSize)


macro(check_functions_exist_and_define func_names)
  foreach(func_name ${func_names})
    string(TOUPPER ${func_name} DEFINE_SUFFIX)
    check_function_exists(${func_name} HAVE_${DEFINE_SUFFIX})
    # message(STATUS "************* ${func_name} ${HAVE_${DEFINE_SUFFIX}} ")
  endforeach()
endmacro(check_functions_exist_and_define)

macro(check_headers_exist_and_define header_names)
  foreach(header_name ${header_names})
    string(TOUPPER ${header_name} DEFINE_SUFFIX)
    string(REPLACE "." "_" DEFINE_SUFFIX ${DEFINE_SUFFIX})
    string(REPLACE "/" "_" DEFINE_SUFFIX ${DEFINE_SUFFIX})
    check_include_file(${header_name} HAVE_${DEFINE_SUFFIX})
    # message(STATUS "************* ${header_name} ${HAVE_${DEFINE_SUFFIX}} ")
  endforeach()
endmacro(check_headers_exist_and_define)


set(FUNCS_TO_CHECK accept memmove memset strdup)
# The quote is needed here to evaluate
check_functions_exist_and_define("${FUNCS_TO_CHECK}")


set(HEADERS_TO_CHECK arpa/nameser.h ctype.h dlfcn.h inttypes.h malloc.h memory.h 
                     netdb.h netinet/in.h netinet/tcp.h resolv.h socket.h
                     stdarg.h stdint.h stdlib.h strings.h string.h sys/malloc.h
                     sys/socket.h sys/stat.h sys/types.h unistd.h)
# The quote is needed here to evaluate
check_headers_exist_and_define("${HEADERS_TO_CHECK}")

check_include_files("stdlib.h;stdarg.h;string.h;float.h" STDC_HEADERS)

# We use this just to get the define. It will check existence as well
# and sets HAVE_INTPTR_T
check_type_size(intptr_t INTPTR_T)

# Generate a configure header
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/antlr3config.h.in.cmake ${CMAKE_CURRENT_SOURCE_DIR}/include/antlr3config.h)



add_library(omantlr3 STATIC)

target_sources(omantlr3 PRIVATE ${ANTLR_SOURCES})
target_include_directories(omantlr3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
