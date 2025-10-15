cmake_minimum_required(VERSION 3.18)

project(primme
  VERSION 3.2.3
  DESCRIPTION "PRIMME: PReconditioned Iterative MultiMethod Eigensolver"
  LANGUAGES C CXX
  HOMEPAGE_URL https://github.com/primme/primme)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

if (POLICY CMP0042)
  # enable MACOSX_RPATH by default
  cmake_policy (SET CMP0042 NEW)
endif ()

# Find our extra modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(PrimmeTests)

find_package(BLAS REQUIRED)
if (NOT TARGET LAPACK::LAPACK)
  find_package(LAPACK REQUIRED)
endif()
primme_test_fortran_underscore()

CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)
if (HAVE_LIB_M)
  set(EXTRA_LIBS ${EXTRA_LIBS} m)
endif()

add_library(primme STATIC)
target_compile_options(primme PRIVATE -fPIC)

target_include_directories(primme
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/include>
)

target_link_libraries(primme PUBLIC LAPACK::LAPACK BLAS::BLAS ${EXTRA_LIBS})

add_subdirectory(include)
add_subdirectory(src)

install(TARGETS primme)

if(BUILD_MY_EXAMPLE)
    add_executable(my_example examples/my_example.c)
    target_include_directories(my_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(my_example PRIVATE primme)
endif()
