cmake_minimum_required(VERSION 3.10)
project(moo)

set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)
set(THIRD_PARTY_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third-party-install/usr/local)

find_package(PkgConfig)

pkg_check_modules(COINHSL coinhsl)

include(ProcessorCount)
ProcessorCount(N)

set(SOURCES
    # base sources for GDOP and similar problems
    src/base/fLGR.cpp
    src/base/trajectory.cpp
    src/base/mesh.cpp
    src/base/linalg.cpp

    # General Dynamic Optimization Problem
    src/nlp/instances/gdop/gdop.cpp
    src/nlp/instances/gdop/problem.cpp
    src/nlp/instances/gdop/strategies.cpp
    src/nlp/instances/gdop/orchestrator.cpp

    # generic NLP
    src/nlp/solvers/nlp_solver_settings.cpp
    src/nlp/nlp_scaling.cpp
    src/nlp/nlp.cpp

    # Ipopt Interface for generic NLP
    src/nlp/solvers/ipopt/adapter.cpp
    src/nlp/solvers/ipopt/solver.cpp

    # C interface
    src/interfaces/c/problem.cpp

    # Sources to call simulations for C interface
    src/simulation/radau/wrapper.cpp
    src/simulation/radau/test.cpp
)

add_library(moo SHARED ${SOURCES})

target_include_directories(moo
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)

###### START FMT

add_subdirectory(third-party/fmt)
target_compile_options(fmt PRIVATE -fPIC)

###### END FMT

###### START DISABLE_LOGGING

option(DISABLE_LOGGING "Disable logging macros" FALSE)

if(DISABLE_LOGGING)
  message(STATUS "DISABLE_LOGGING set to TRUE")
  target_compile_definitions(moo PRIVATE DISABLE_LOGGING)
else()
  message(STATUS "DISABLE_LOGGING set to FALSE")
endif()

###### END DISABLE_LOGGING

###### START MUMPS + Ipopt

# TODO: enable build w/o mumps or ipopt - alternative NLP solvers
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps/MUMPS)
  message(STATUS "Downloading mumps sources...")
  execute_process(COMMAND ./get.Mumps WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps OUTPUT_QUIET)
endif()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps)
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/config.status
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps/configure
        --enable-shared=no
        --enable-static=yes
  COMMAND make clean
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps)
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs/libcoinmumps.a
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/config.status
  COMMAND make -j
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps)
add_custom_target(coinmumps-build DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs/libcoinmumps.a)
add_library(moo::coinmumps STATIC IMPORTED)
add_dependencies(moo::coinmumps coinmumps-build)
set_property(TARGET moo::coinmumps PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs/libcoinmumps.a)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt)
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/config.status
  COMMAND ADD_CFLAGS=-Ofast ADD_CXXFLAGS=-Ofast ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt/configure
      --without-pardiso
      --with-mumps
      --with-mumps-lflags='-L${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs'
      --with-mumps-cflags='-DCOIN_USE_MUMPS_MPI_H -I${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps -I${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps/MUMPS/include -I${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps/MUMPS/libseq'
      --enable-shared=no
      --enable-static=yes
  COMMAND make clean
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt)
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/src/.libs/libipopt.a
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/config.status
  COMMAND make -j
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt)
add_custom_target(ipopt-build DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/src/.libs/libipopt.a)
add_library(moo::ipopt STATIC IMPORTED)
add_dependencies(moo::ipopt ipopt-build)
set_property(TARGET moo::ipopt PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/src/.libs/libipopt.a)
target_include_directories(moo::ipopt INTERFACE 
  ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt/src/Interfaces
  ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt/src/Common
  ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt/src/Algorithm
  ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt/src/LinAlg
)

###### END MUMPS + Ipopt

###### START RADAU5

option(WITH_RADAU "Enable RADAU build and linking" ON)

if(WITH_RADAU)
    add_subdirectory(third-party/RADAU)
endif()

##### END RADAU5

target_include_directories(moo PUBLIC include)
target_include_directories(moo PUBLIC src)

target_include_directories(moo PRIVATE ${THIRD_PARTY_INSTALL_PREFIX}/include/coin-or)

target_link_directories(moo PRIVATE ${THIRD_PARTY_INSTALL_PREFIX}/lib)

#add_compile_options(-fsanitize=undefined -fsanitize=address -fstack-protector-strong)
#add_link_options(-fsanitize=undefined -fsanitize=address -fstack-protector-strong)

# TODO: Make it a clean build with -Wextra

# force release build if included in 3rdParty tool - explicitly set MOO_BUILD_TYPE to debug for debug build!
if(NOT DEFINED MOO_BUILD_TYPE)
  set(MOO_BUILD_TYPE "RelWithDebInfo")
  message(STATUS "MOO_BUILD_TYPE not defined: defaulting to RelWithDebInfo")
endif()

if ("${MOO_BUILD_TYPE}" STREQUAL "Debug")
    message(STATUS "Configuration: Debug")
    target_compile_options(moo PRIVATE -Wall -fno-rtti -fno-var-tracking-assignments -Wnon-virtual-dtor -pedantic -O0 -g)
    message(STATUS "Compile Options: -Wall -fno-rtti -fno-var-tracking-assignments -Wnon-virtual-dtor -pedantic -O0 -g")
else ()
    message(STATUS "Configuration: RelWithDebInfo / Release")
    target_compile_options(moo PRIVATE -Wall -fno-rtti -fno-var-tracking-assignments -Wnon-virtual-dtor -pedantic -O3 -march=native -funroll-loops -g)
    message(STATUS "Compile Options: -Wall -fno-rtti -fno-var-tracking-assignments -Wnon-virtual-dtor -pedantic -O3 -march=native -funroll-loops -g")
endif()

target_link_libraries(moo PUBLIC fmt)
target_link_libraries(moo PRIVATE moo::ipopt moo::coinmumps metis)
target_link_libraries(moo PRIVATE lapack blas dl gfortran gomp)

if(WITH_RADAU)
    target_link_libraries(moo PRIVATE radau)
endif()

# hide symbols
target_compile_options(moo PRIVATE -Wl,--exclude-libs,ALL )
target_link_options(moo PRIVATE -Wl,--exclude-libs,ALL)

if (COINHSL_FOUND)
    message(STATUS "Coin-HSL found")
    target_link_libraries(moo PRIVATE coinhsl)
endif ()

install(TARGETS moo
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
