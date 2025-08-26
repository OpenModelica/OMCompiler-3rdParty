cmake_minimum_required(VERSION 3.10)

option(WITH_RADAU "Build and link with 'RADAU' library" TRUE)

if(WITH_RADAU)
    message(STATUS "Library 'RADAU' will be built from source and linked (WITH_RADAU=TRUE)")

    project(radau)

    enable_language(Fortran)

    find_package(LAPACK REQUIRED)

    add_library(radau STATIC
        radau.f
        dc_lapack.f
    )

    target_compile_options(radau PUBLIC -O3 -fPIC ${MOO_NATIVE_BUILD_FLAGS})

    target_link_libraries(radau PUBLIC ${LAPACK_LIBRARIES} gfortran)

    target_include_directories(radau PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

else()
  message(STATUS "Library 'WITH_RADAU' is disabled (WITH_RADAU=FALSE)")
endif()

