cmake_minimum_required(VERSION 3.4...3.14)
project(zlib)

# zlib
# We have decided to use zlib from here. We could have used the system zlib. However,
# modelica annotations request for "zlib" while the system zlib is OFTEN (but not always) called "libz"
# which means it should be used as "z". We can modify the annotations to use "z" but then
# it will be the same issue on systems that call it "zlib". So we need to find a solution.
# Originally I was creating a sym link to the system zlib in our lib directories. However,
# that might be confusing for others. So it might be better to explicitly
# build it and use it from here. The one advantage of this is that we can compile it with -fpic so
# that we can link it into our static FMUs.
file(GLOB libzlib_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
add_library(zlib STATIC ${libzlib_SOURCES})
add_library(omc::3rd::zlib ALIAS zlib)


target_include_directories(zlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

## Check if we have unistd.h and define accordingly
include(CheckIncludeFile)
check_include_file(unistd.h HAVE_UNISTD_H)
if(HAVE_UNISTD_H)
  target_compile_definitions(zlib PRIVATE -DHAVE_UNISTD_H=1)
endif()

install(TARGETS zlib
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/zlib.h
              ${CMAKE_CURRENT_SOURCE_DIR}/zconf.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
