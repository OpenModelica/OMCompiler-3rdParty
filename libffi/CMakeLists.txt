# This CMake file is originally taken from https://github.com/am11/libffi/tree/feature/cmake-build-configs with love.
# It has been heavily modified to adapt to more recent and in some places correct cmake conventions.
# It is intended to be used with the OpenModelica CMake build system.

cmake_minimum_required(VERSION 3.14)
project(libffi C)

option(LIBFFI_BUILD_STATIC_LIBS "Build libffi as static lib." ON)
option(LIBFFI_BUILD_SHARED_LIBS "Build libffi as shared lib." OFF)

if(LIBFFI_BUILD_STATIC_LIBS AND LIBFFI_BUILD_SHARED_LIBS)
    message(FATAL_ERROR "You have asked for both shared and static versions of Libffi. This is not recommended. Please pick one.")
endif()


if(LIBFFI_BUILD_STATIC_LIBS)
    add_library(ffi STATIC)
elseif(LIBFFI_BUILD_SHARED_LIBS)
    add_library(ffi SHARED)
else()
    message(FATAL_ERROR "You have disabled both shared and static versions of Libffi. Please pick one.")
endif()

set(SOURCES_LIST
    src/closures.c
    src/java_raw_api.c
    src/prep_cif.c
    src/raw_api.c
    src/types.c)

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
    list(APPEND SOURCES_LIST src/debug.c)
    target_compile_definitions(ffi PRIVATE -DFFI_DEBUG)
endif()

include(.cmake/configure_platform.cmake)
include(.cmake/configure_options.cmake)

target_sources(ffi PRIVATE ${SOURCES_LIST})

set_property(TARGET ffi PROPERTY POSITION_INDEPENDENT_CODE 1)

target_include_directories(ffi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# Some configured files are generated in the build directory.
target_include_directories(ffi PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)

target_compile_definitions(ffi PRIVATE -DFFI_BUILDING)



if(LIBFFI_BUILD_SHARED_LIBS)
    target_sources(ffi PRIVATE ${OBJECTS_LIST})

    if(MSVC)
        target_compile_definitions(ffi PRIVATE -DFFI_BUILDING_DLL)
        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    elseif(CMAKE_SYSTEM_NAME MATCHES Darwin)
        set(CMAKE_MACOSX_RPATH 1)
    endif()
endif()

# add_library(ffi_shared SHARED $<TARGET_OBJECTS:objlib> ${OBJECTS_LIST})

file(COPY src/${TARGETDIR}/ffitarget.h DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include)

install(TARGETS ffi
    EXPORT ${PROJECT_NAME}Targets)

# install(EXPORT ${PROJECT_NAME}Targets
#     DESTINATION share/${PROJECT_NAME})

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/include/ffi.h
    ${CMAKE_CURRENT_BINARY_DIR}/include/ffitarget.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ffi)
