cmake_minimum_required(VERSION 3.14)

if(MOO_WITH_RADAU)
    message(STATUS "Library RADAU will be built from source and linked (MOO_WITH_RADAU=TRUE)")

    project(radau)

    enable_language(Fortran)

    add_library(radau STATIC
        radau.f
        dc_lapack.f
    )

    target_compile_options(radau PRIVATE -O3 -fPIC)

    target_link_libraries(radau PUBLIC ${MOO_LAPACK_LIB} ${MOO_GFORTRAN_LIB})

    target_include_directories(radau PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

else()
  message(STATUS "Library 'MOO_WITH_RADAU' is disabled (MOO_WITH_RADAU=FALSE)")
endif()

