cmake_minimum_required(VERSION 3.10)
project(moo LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

##### START setup

message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")

if(WIN32)
    set(LIB_EXT .dll)
elseif(APPLE)
    set(LIB_EXT .dylib)
else()
    set(LIB_EXT .so)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(IS_GCC TRUE)
else()
    set(IS_GCC FALSE)
endif()

##### END setup

##### START make

option(MOO_MAKEPROGRAM "Make program used to build third party dependencies of MOO library")

# A make program is necessary for the legacy configure builds of Ipopt and MUMPS
if(${MOO_MAKEPROGRAM})
    set(MAKEPROGRAM "${MOO_MAKEPROGRAM}")
else()
    if($ENV{MAKE})
        set(MAKEPROGRAM "$ENV{MAKE}")
    else()
        # Now we might be in trouble
        find_program(MAKEPROGRAM "make")
        if("${MAKEPROGRAM}" STREQUAL "MAKEPROGRAM-NOTFOUND") # last ditch
            set(MAKEPROGRAM "$$(MAKE)")
            message(WARNING "Could not find a make program. Falling back to $(MAKE) enviroment variable. Please set MOO_MAKEPROGRAM to specify an appropriate make program.")
            if("${CMAKE_MAKE_PROGRAM}" MATCHES "ninja$")
                message(SEND_ERROR "Ninja build does not support reading $(MAKE) environment variable. This will probably fail.")
            endif()
        endif()
    endif()
endif()
message(STATUS "Make program used to build third party dependencies of MOO library: ${MAKEPROGRAM}")

##### END make

include(ExternalProject)
set(THIRD_PARTY_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third-party-install/usr/local)

find_package(PkgConfig)

###### START DEPENDENCIES

pkg_check_modules(COINHSL coinhsl)

set(MOO_LAPACK_LIB "" CACHE STRING "Override path to LAPACK library")
set(MOO_METIS_LIB "" CACHE STRING "Override path to METIS library")
set(MOO_GFORTRAN_LIB "" CACHE STRING "Override path to gfortran library")

set(MOO_EXTRA_LIBS "")

# LAPACK
if(MOO_LAPACK_LIB)
    message(STATUS "Using user-provided LAPACK: ${MOO_LAPACK_LIB}")
    list(APPEND MOO_EXTRA_LIBS ${MOO_LAPACK_LIB})
else()
    message(STATUS "MOO_LAPACK_LIB not defined: searching for LAPACK")
    find_package(LAPACK REQUIRED)

    if(LAPACK_FOUND)
        set(MOO_LAPACK_LIB ${LAPACK_LIBRARIES})
        message(STATUS "Found LAPACK: ${MOO_LAPACK_LIB}")
        list(APPEND MOO_EXTRA_LIBS ${MOO_LAPACK_LIB})
    else()
        message(FATAL_ERROR "LAPACK library not found")
    endif()
endif()

# METIS
if(MOO_METIS_LIB)
    message(STATUS "Using user-provided METIS: ${MOO_METIS_LIB}")
    list(APPEND MOO_EXTRA_LIBS ${MOO_METIS_LIB})
else()
    message(STATUS "MOO_METIS_LIB not defined: searching for METIS")
    find_library(MOO_METIS_LIB NAMES metis)
    if(MOO_METIS_LIB)
        message(STATUS "Found METIS: ${MOO_METIS_LIB}")
        list(APPEND MOO_EXTRA_LIBS ${MOO_METIS_LIB})
    else()
        message(FATAL_ERROR "METIS library not found")
    endif()
endif()

# gfortran runtime
if(MOO_GFORTRAN_LIB)
    message(STATUS "Using user-provided gfortran: ${MOO_GFORTRAN_LIB}")
    list(APPEND MOO_EXTRA_LIBS ${MOO_GFORTRAN_LIB})
else()
    message(STATUS "MOO_GFORTRAN_LIB not defined: searching for gfortran")
    find_library(MOO_GFORTRAN_LIB NAMES gfortran)

    if(NOT MOO_GFORTRAN_LIB)
        message(STATUS "gfortran library not found in default paths. Attempting to find via the compiler.")
        find_program(GFORTRAN_EXECUTABLE NAMES gfortran)

        if(GFORTRAN_EXECUTABLE)
            message(STATUS "Found gfortran executable: ${GFORTRAN_EXECUTABLE}")

            execute_process(
                COMMAND "${GFORTRAN_EXECUTABLE}" -print-file-name=libgfortran${LIB_EXT}
                OUTPUT_VARIABLE GFORTRAN_LIB_PATH
                OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE EXECUTE_RESULT
            )

            if(GFORTRAN_LIB_PATH)
                message(STATUS "Found gfortran library path using the compiler: ${GFORTRAN_LIB_PATH}")
                set(MOO_GFORTRAN_LIB "${GFORTRAN_LIB_PATH}")
            else()
                message(STATUS "Failed to get gfortran library path from the compiler.")
            endif()
        else()
            message(STATUS "gfortran executable not found.")
        endif()
    else()
        message(STATUS "Found gfortran: ${MOO_GFORTRAN_LIB}")
        list(APPEND MOO_EXTRA_LIBS ${MOO_GFORTRAN_LIB})
    endif()

    if(NOT MOO_GFORTRAN_LIB)
        message(WARNING "gfortran library not found: linking may still work if system finds it automatically")
    endif()
endif()

###### END DEPENDENCIES


###### START BUILD + COMPILER FLAGS

# Attention:
# force release build if included in 3rdParty tool - explicitly set MOO_BUILD_TYPE to debug for debug build!
if(NOT DEFINED MOO_BUILD_TYPE)
  set(MOO_BUILD_TYPE "RelWithDebInfo")
  message(STATUS "MOO_BUILD_TYPE not defined: defaulting to RelWithDebInfo")
endif()

if(MOO_NATIVE_BUILD)
    set(MOO_NATIVE_BUILD_FLAGS
        -march=native
        -mtune=native
        -funroll-loops
        -ftree-vectorize
        -fprefetch-loop-arrays
    )
    string(JOIN " " MOO_NATIVE_BUILD_FLAGS_STRING ${MOO_NATIVE_BUILD_FLAGS})
    message(STATUS "Adding native CPU optimizations: ${MOO_NATIVE_BUILD_FLAGS_STRING}")

else()
    set(MOO_NATIVE_BUILD_FLAGS_STRING "")
    message(STATUS "MOO_NATIVE_BUILD not defined: adding no additional compiler flags")
endif()

###### END BUILD + COMPILER FLAGS


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

###### START MUMPS + Ipopt

# TODO: enable build w/o mumps or ipopt - alternative NLP solvers
# TODO: remove Ipopt examples / tutorial as they bloat the build

if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps/MUMPS)
  message(STATUS "Downloading mumps sources...")
  execute_process(COMMAND ./get.Mumps WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps OUTPUT_QUIET)
endif()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/config.status
    COMMAND CC=${CMAKE_C_COMPILER}
            CXX=${CMAKE_CXX_COMPILER}
            FC=${CMAKE_Fortran_COMPILER}
            ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps/configure
            ADD_CFLAGS=-O3\ ${MOO_NATIVE_BUILD_FLAGS_STRING}
            ADD_FCFLAGS=-O3\ ${MOO_NATIVE_BUILD_FLAGS_STRING}
            ADD_CXXFLAGS=-O3\ ${MOO_NATIVE_BUILD_FLAGS_STRING}
            --with-lapack="${MOO_LAPACK_LIB}"
            --enable-shared=no
            --enable-static=yes
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps
    COMMENT "Configuring MUMPS library..."
)

add_custom_target(configure_mumps
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/config.status
    COMMENT "Ensuring MUMPS is configured."
)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs/libcoinmumps.a
    DEPENDS configure_mumps
    COMMAND ${MAKEPROGRAM}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps
    COMMENT "Building MUMPS library..."
)

add_custom_target(coinmumps-build DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs/libcoinmumps.a)
add_library(moo::coinmumps STATIC IMPORTED)
add_dependencies(moo::coinmumps coinmumps-build)
set_property(TARGET moo::coinmumps PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs/libcoinmumps.a)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/config.status
    COMMAND CC=${CMAKE_C_COMPILER}
            CXX=${CMAKE_CXX_COMPILER}
            FC=${CMAKE_Fortran_COMPILER}
            ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt/configure
            ADD_CFLAGS=-O3\ ${MOO_NATIVE_BUILD_FLAGS_STRING}
            ADD_FCFLAGS=-O3\ ${MOO_NATIVE_BUILD_FLAGS_STRING}
            ADD_CXXFLAGS=-O3\ ${MOO_NATIVE_BUILD_FLAGS_STRING}
            --with-lapack="${MOO_LAPACK_LIB}"
            --without-pardiso
            --with-mumps
            --with-mumps-lflags='-L${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs'
            --with-mumps-cflags='-DCOIN_USE_MUMPS_MPI_H -I${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps -I${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps/MUMPS/include -I${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps/MUMPS/libseq'
            --enable-shared=no
            --enable-static=yes
    COMMAND ${MAKEPROGRAM} clean
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt
    COMMENT "Configuring IPOPT library..."
)

add_custom_target(configure_ipopt
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/config.status
    COMMENT "Ensuring IPOPT is configured."
)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/src/.libs/libipopt.a
    DEPENDS configure_ipopt
    COMMAND ${MAKEPROGRAM}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt
    COMMENT "Building IPOPT library..."
)

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

set(MOO_COMPILE_OPTIONS -Wall -fno-rtti -Wnon-virtual-dtor -pedantic)

if(IS_GCC)
    list(APPEND MOO_COMPILE_OPTIONS -fno-var-tracking-assignments)
endif()

if("${MOO_BUILD_TYPE}" STREQUAL "Debug")
    message(STATUS "Configuration: Debug")
    list(APPEND MOO_COMPILE_OPTIONS -O0 -g)
else()
    message(STATUS "Configuration: RelWithDebInfo / Release")
    list(APPEND MOO_COMPILE_OPTIONS -O3 -g ${MOO_NATIVE_BUILD_FLAGS})
endif()

target_compile_options(moo PRIVATE ${MOO_COMPILE_OPTIONS})
message(STATUS "Compile Options: ${MOO_COMPILE_OPTIONS}")

target_link_libraries(moo PUBLIC fmt)
target_link_libraries(moo PRIVATE moo::ipopt moo::coinmumps)
target_link_libraries(moo PRIVATE ${MOO_EXTRA_LIBS})

if(WITH_RADAU)
    target_link_libraries(moo PRIVATE radau)
endif()

if (COINHSL_FOUND)
    target_link_libraries(moo PRIVATE coinhsl)
endif ()

# hide symbols
if(WIN32)
    target_compile_definitions(moo PRIVATE MOO_DLL_EXPORT)
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    target_link_libraries(moo PRIVATE dl)
    target_link_directories(moo PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs/)
    target_link_directories(moo PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/src/.libs/)
    target_link_options(moo PRIVATE -Wl,-hidden-lcoinmumps -Wl,-hidden-lipopt)
    target_compile_options(moo PRIVATE -fvisibility=hidden )
else()
    target_link_libraries(moo PRIVATE dl)
    target_link_options(moo PRIVATE -Wl,--exclude-libs,ALL)
    target_compile_options(moo PRIVATE -fvisibility=hidden)
endif()

install(TARGETS moo
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
