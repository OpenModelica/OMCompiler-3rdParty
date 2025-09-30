cmake_minimum_required(VERSION 3.14)

if(NOT MOO_WITH_GDOPT)
    message(FATAL_ERROR "-DMOO_WITH_GDOPT=ON is required to build tests (-DMOO_TESTS=ON).")
endif()

if(NOT MOO_WITH_RADAU)
    message(FATAL_ERROR "-DMOO_WITH_RADAU=ON is required to build tests (-DMOO_TESTS=ON).")
endif()

add_subdirectory(gdopt)
