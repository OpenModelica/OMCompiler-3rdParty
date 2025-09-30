cmake_minimum_required(VERSION 3.14)
project(Ipopt)

# FUNNY_LAPACK_FINT

#undef HAVE_MPI_INITIALIZED - Define to 1 if you have the `MPI_Initialized' function.

option(IPOPT_HAS_PARDISO                "Enable Pardiso solver" OFF)
option(IPOPT_HAS_PARDISO_MKL            "Enable if you are using Pardiso from MKL" OFF)
option(IPOPT_HAS_PARDISO_PARALLEL       "Enable if you are using the parallel version of Pardiso" OFF)
option(IPOPT_HAS_HSL                    "Enable HSL interface"  OFF)
option(IPOPT_HAS_WSMP                   "Enable WSMP solver"    OFF)
option(IPOPT_HAS_SPRAL                  "Enable SPRAL solver"   OFF)
option(IPOPT_HAS_MUMPS                  "Enable Mumps solver"   ON)
option(IPOPT_SINGLE                     "Build with single precision"     OFF)
option(IPOPT_INT64                      "Build with 64 bit integer type"  OFF)
option(IPOPT_BUILD_EXAMPLES             "Enable the building of examples" OFF)
option(IPOPT_HAS_LINEARSOLVERLOADER     "Build the dynamic linear solver loader"  OFF)
option(IPOPT_ENABLE_PARDISOSOLVERLOADER "Build the dynamic pardiso solver loader" OFF)
option(IPOPT_ENABLE_INEXACT             "Build the inexact solver" OFF)

set(IPOPT_CHECKLEVEL "0" CACHE STRING "The debug sanity check level of Ipopt (0 if no test)")
set(IPOPT_VERBOSITY  "0" CACHE STRING "The debug verbosity level of Ipopt (0 if no output)")

# options from ThirdParty/CMakeLists.txt
option(DOWNLOAD_MUMPS               "Enable the download / compilation of Mumps" ON)
option(DOWNLOAD_LAPACK              "Enable the download / compilation of Blas / Lapack" OFF)
option(USE_SYSTEM_LAPACK            "Enable the use of the system Lapack" ON)
option(USE_PROCESSOR_EXTENSIONS     "Use sse / mmx / avx extensions during compilation" OFF)

# MUMPS build options
option(MUMPS_USE_LIBSEQ     "Use the MUMPS sequential MPI stub" ON)
option(MUMPS_USE_METIS      "Use the Metis library"             ON)

# TODO: add OpenMP option
option(MUMPS_USE_OPENMP     "Use the OpenMP library"            OFF)

# TODO: add MPI option
option(MUMPS_USE_MPI        "Use the MPI library"               OFF)

set(METIS_LIB_PATH "" CACHE STRING "Path to libmetis (file) or directory containing libmetis")
set(METIS_INC_PATH "" CACHE PATH   "Path to METIS include directory (contains metis.h)")
set(EXTERNAL_METIS_TARGET "" CACHE STRING "Existing METIS target to use: add_library(EXTERNAL_METIS_TARGET ALIAS MY_METIS_TARGET)")

# Compilation options
option(USE_FAST_CODE            "Use intensive optimization flags" ON)
option(ENABLE_SHARED_LIBRARIES  "Build libraries as shared libraries" OFF)
option(COMPILE_LTO              "Activate the whole program optimization" OFF)  # GCC macOS does not support this
option(HAS_MKL                  "Use Intel MKL library (requires Intel compiler)" OFF)

# Set paths of source
set(Ipopt_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Ipopt/sources/)

mark_as_advanced(ENABLE_SHARED_LIBRARIES
        COMPILE_LTO
        HAS_MKL
        USE_FAST_CODE
        IPOPT_HAS_PARDISO_MKL
        IPOPT_HAS_PARDISO_PARALLEL
        IPOPT_CHECKLEVEL
        IPOPT_VERBOSITY
        USE_PROCESSOR_EXTENSIONS)

# Set paths for binary and library generation inside the build directory:
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(GetAcInitVersion)
include(GNUInstallDirs)
include(cmake/utils.cmake)

# TODO: remove this awful solution by explicitly setting the flags for targets instead of globally!!
safe_tmp(_GLOBL_TMP_CFLAGS       CMAKE_C_FLAGS)
safe_tmp(_GLOBL_TMP_CXXFLAGS     CMAKE_CXX_FLAGS)
safe_tmp(_GLOBL_TMP_FortranFLAGS CMAKE_Fortran_FLAGS)

Enable_Testing ()

#------------------------------------------------------------
# Check options
#------------------------------------------------------------

if (IPOPT_ENABLE_INEXACT AND NOT (IPOPT_HAS_PARDISO OR IPOPT_HAS_PARDISO_MKL OR IPOPT_HAS_PARDISO_PARALLEL))
    message(FATAL_ERROR "Error: Inexact solver is only available through MKL. Please activate the MKL")
endif ()

#------------------------------------------------------------
# Detect 64 bits
#------------------------------------------------------------

if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(HAVE_64_BIT 0)
else ()
    set(HAVE_64_BIT 1)
endif ()

# Various definitions

# Name of package
set(PACKAGE           "Ipopt")
# Define to the address where bug reports for this package should be sent.
set(PACKAGE_BUGREPORT "bugs@coin-or.org")
# Define to the full name of this package.
set(PACKAGE_NAME      "Ipopt")
# Define to the full name and version of this package.
set(PACKAGE_STRING    "Ipopt")
# Define to the one symbol short name of this package.
set(PACKAGE_TARNAME   "ipopt")

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/install" CACHE PATH "Ipopt install prefix" FORCE)
endif ()

include(CheckCCompilerFlag)
if (COMPILE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported()
endif ()

if (USE_PROCESSOR_EXTENSIONS)
    # Check for SSE* and AVX*
    find_package(SSE)
    if (MMX_FOUND OR SSE2_FOUND OR SSE3_FOUND OR SSSE3_FOUND OR
        SSE4_1_FOUND OR SSE4_2_FOUND OR AVX_FOUND OR AVX2_FOUND)
        set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} ${SSE_COMPILER_FLAGS}")
        set(CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS} ${SSE_COMPILER_FLAGS}")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${SSE_COMPILER_FLAGS}")
    endif ()
endif ()

add_library(METIS_TARGET INTERFACE)

if (MUMPS_USE_METIS)
    if(TARGET EXTERNAL_METIS_TARGET)
        message(STATUS "Using external METIS target.")
        target_link_libraries(METIS_TARGET INTERFACE EXTERNAL_METIS_TARGET)
    elseif(METIS_LIB_PATH AND METIS_INC_PATH)
        add_library(_METIS_USER_IMPORTED UNKNOWN IMPORTED)
        set_target_properties(_METIS_USER_IMPORTED PROPERTIES
            IMPORTED_LOCATION "${METIS_LIB_PATH}"
            INTERFACE_INCLUDE_DIRECTORIES "${METIS_INC_PATH}"
        )
        message(STATUS "Using user-provided METIS: ${METIS_LIB_PATH}")
        target_link_libraries(METIS_TARGET INTERFACE _METIS_USER_IMPORTED)
    else()
        find_package(METIS)
        if(NOT METIS_FOUND)
            message(FATAL_ERROR "METIS is required for MUMPS but not found. "
                                "Provide -DMETIS_LIB_PATH=/path/to/libmetis.so "
                                "and -DMETIS_INC_PATH=/path/to/include, "
                                "or disable METIS with -DMUMPS_USE_METIS=OFF")
        else()
            target_link_libraries(METIS_TARGET INTERFACE METIS::METIS)
        endif()
    endif()
endif()

# Check for MKL
if (HAS_MKL)
    find_package(MKL)

    if (MKL_FOUND)
        message(STATUS "MKL library found")
    else ()
        message(STATUS "MKL library not found")
    endif ()

    # Copy libiomp5md.dll in the build directory
    if (WIN32)
        if (HAVE_64_BIT)
            set(MKL_DLL_DIR ${MKL_ROOT}/bin/intel64)
        else ()
            set(MKL_DLL_DIR ${MKL_ROOT}/bin/ia32)
        endif ()

        execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/bin/
                COMMAND ${CMAKE_COMMAND} -E copy ${MKL_DLL_DIR}/libiomp5md.dll ${CMAKE_CURRENT_BINARY_DIR}/bin/
                COMMAND ${CMAKE_COMMAND} -E echo "Copying ${MKL_DLL_DIR}/libiomp5md.dll into ${CMAKE_CURRENT_BINARY_DIR}/bin/")
    endif ()

    set(COIN_MKL_LIBS "${MKL_LIBRARIES}")
    if (WIN32)
        set(COIN_MKL_LIBS ${COIN_MKL_LIBS} mkl_intel_thread libiomp5md)
    else ()
        set(COIN_MKL_LIBS ${COIN_MKL_LIBS} mkl_gnu_thread gomp dl)
    endif ()

    include_directories(${MKL_INCLUDE_DIRS})

    if (HAVE_64_BIT)
        link_directories(${MKLROOT_PATH}/mkl/lib/intel64)
    else ()
        link_directories(${MKLROOT_PATH}/mkl/lib/ia32)
    endif ()

    set(COIN_HAS_LAPACK ON CACHE BOOL "Use Intel MKL library (requires Intel compiler)")
endif ()

#-----------------------------------------------------------------------------
# Dependencies
#-----------------------------------------------------------------------------

if (USE_SYSTEM_LAPACK AND NOT DOWNLOAD_LAPACK)
    find_package(LAPACK REQUIRED)
endif ()


include(CheckIncludeFileCXX)
include(CheckIncludeFile)
include(VA_COPY)

check_include_file("assert.h"    HAVE_ASSERT_H)
check_include_file("bzlib.h"     HAVE_BZLIB_H)
check_include_file("ctype.h"     HAVE_CTYPE_H)
check_include_file("dlfcn.h"     HAVE_DLFCN_H)
check_include_file("endian.h"    HAVE_ENDIAN_H)
check_include_file("float.h"     HAVE_FLOAT_H)
check_include_file("ieeefp.h"    HAVE_IEEEFP_H)
check_include_file("inttypes.h"  HAVE_INTTYPES_H)
check_include_file("math.h"      HAVE_MATH_H)
check_include_file("memory.h"    HAVE_MEMORY_H)
check_include_file("stdint.h"    HAVE_STDINT_H)
check_include_file("stdlib.h"    HAVE_STDLIB_H)
check_include_file("stdio.h"     HAVE_STDIO_H)
check_include_file("stdarg.h"    HAVE_STDARG_H)
check_include_file("stddef.h"    HAVE_STDDEF_H)
check_include_file("strings.h"   HAVE_STRINGS_H)
check_include_file("string.h"    HAVE_STRING_H)
check_include_file("sys/stat.h"  HAVE_SYS_STAT_H)
check_include_file("sys/types.h" HAVE_SYS_TYPES_H)
check_include_file("time.h"      HAVE_TIME_H)
check_include_file("unistd.h"    HAVE_UNISTD_H)
check_include_file("windows.h"   HAVE_WINDOWS_H)
check_include_file("zlib.h"      HAVE_ZLIB_H)
check_include_file_cxx("cctype"    HAVE_CCTYPE)
check_include_file_cxx("cmath"     HAVE_CMATH)
check_include_file_cxx("cieeefp"   HAVE_CIEEEFP)
check_include_file_cxx("cfloat"    HAVE_CFLOAT)
check_include_file_cxx("cinttypes" HAVE_CINTTYPES)
check_include_file_cxx("cassert"   HAVE_CASSERT)
check_include_file_cxx("cstdio"    HAVE_CSTDIO)
check_include_file_cxx("cstdlib"   HAVE_CSTDLIB)
check_include_file_cxx("cstdarg"   HAVE_CSTDARG)
check_include_file_cxx("cstddef"   HAVE_CSTDDEF)
check_include_file_cxx("cstring"   HAVE_CSTRING)
check_include_file_cxx("ctime"     HAVE_CTIME)

string(SUBSTRING ${CMAKE_SHARED_LIBRARY_SUFFIX} 1 -1 SHAREDLIBEXT)

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag(-Qunused-arguments HAVE_QUNUSED_ARGUMENTS)
if (HAVE_QUNUSED_ARGUMENTS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")
endif ()

if (UNIX)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT APPLE)
        check_cxx_compiler_flag(-Wno-unused-local-typedefs GCC_HAS_TYPEDEFS)
        if (GCC_HAS_TYPEDEFS)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
        endif ()
    endif ()

    if (NOT APPLE)
        check_cxx_compiler_flag(-Wno-narrowing GCC_HAS_NARROWING)
        if (GCC_HAS_NARROWING)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")
            set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Wno-narrowing")
        endif ()
    endif ()
endif ()

include(CheckTypeSize)

check_type_size("int *" SIZEOF_INT_P)

include(CheckSymbolExists)
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)
include(CheckFunctionExists)

set(CMAKE_REQUIRED_LIBRARIES m)
check_cxx_symbol_exists(std::isfinite "cmath" IPOPT_C_FINITE_R)
if (IPOPT_C_FINITE_R)
    set(IPOPT_C_FINITE std::isfinite)
endif ()
if (NOT IPOPT_C_FINITE_R)
    set(CMAKE_REQUIRED_LIBRARIES m)
    check_symbol_exists(isfinite "math.h" IPOPT_C_FINITE_R)
    if (IPOPT_C_FINITE_R)
        set(IPOPT_C_FINITE isfinite)
    endif ()
endif ()
if (NOT IPOPT_C_FINITE_R)
    set(CMAKE_REQUIRED_LIBRARIES m)
    check_symbol_exists(_finite "math.h" IPOPT_C_FINITE_R)
    if (IPOPT_C_FINITE_R)
        set(IPOPT_C_FINITE _finite)
    endif ()
endif ()
if (NOT IPOPT_C_FINITE_R)
    set(CMAKE_REQUIRED_LIBRARIES m)
    check_symbol_exists(finite "math.h" IPOPT_C_FINITE_R)
    if (IPOPT_C_FINITE_R)
        set(IPOPT_C_FINITE finite)
    endif ()
endif ()

set(CMAKE_REQUIRED_LIBRARIES m)
check_cxx_symbol_exists(std::isnan "cmath" COIN_C_ISNAN_R)
if (COIN_C_ISNAN_R)
    set(COIN_C_ISNAN std::isnan)
endif ()
if (NOT COIN_C_ISNAN_R)
    set(CMAKE_REQUIRED_LIBRARIES m)
    check_symbol_exists(isnan "math.h" COIN_C_ISNAN_R)
    if (COIN_C_ISNAN_R)
        set(COIN_C_NAN isnan)
    endif ()
endif ()
if (NOT COIN_C_ISNAN_R)
    set(CMAKE_REQUIRED_LIBRARIES m)
    check_symbol_exists(_isnan "math.h" COIN_C_ISNAN_R)
    if (COIN_C_NAN_R)
        set(COIN_C_NAN _isnan)
    endif ()
endif ()
if (NOT COIN_C_ISNAN_R)
    set(CMAKE_REQUIRED_LIBRARIES m)
    check_symbol_exists(isnand "ieeefp.h" COIN_C_ISNAN_R)
    if (COIN_C_NAN_R)
        set(COIN_C_NAN isnand)
    endif ()
endif ()

check_symbol_exists(clock_gettime time.h     HAVE_CLOCK_GETTIME)
check_symbol_exists(gettimeofday  sys/time.h HAVE_GETTIMEOFDAY)

check_function_exists(dran48 IPOPT_HAS_DRAND48)
check_function_exists(rand IPOPT_HAS_RAND)

check_include_file_cxx(cstdlib IPOPT_HAS_STD__RAND)
check_symbol_exists(va_copy "stdarg.h" IPOPT_HAS_VA_COPY)
check_function_exists(fopen_s IPOPT_HAS_FOPEN_S)
check_function_exists(getenv_s IPOPT_HAS_GETENV_S)
check_symbol_exists(feenableexcept "cfenv" IPOPT_HAS_FEENABLEEXCEPT)
check_symbol_exists(sigaction "signal.h" IPOPT_HAS_SIGACTION)

if (IPOPT_ENABLE_INEXACT)
    set(BUILD_INEXACT 1)
endif ()

if (COIN_HAS_LAPACK OR DOWNLOAD_LAPACK OR USE_SYSTEM_LAPACK)
    set(IPOPT_HAS_LAPACK 1)
endif ()

if (WIN32)
    # snprintf not correctly detected under Visual Studio.
    # Hack: we just activate snprintf under Windows.
    # TO BE FIXED
    set(HAVE_STDIO_H    1)
    set(HAVE__SNPRINTF  1)
    set(HAVE_VSNPRINTF  1)
    set(HAVE__VSNPRINTF 1)
    set(IPOPT_HAS_VA_COPY 1)
endif ()

#-----------------------------------------------------------------------------
# Manage compilation options
#-----------------------------------------------------------------------------

if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Wno-deprecated")

    if (NOT ENABLE_SHARED_LIBRARIES)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fPIC")
    endif ()
endif ()

#-----------------------------------------------------------------------------
# Manage specific VS flags
#-----------------------------------------------------------------------------

if (MSVC)
    # Avoid Warning C4530 by using the flag /EHsc
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc ")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} /EHsc ")
endif ()

#-----------------------------------------------------------------------------
# Manage threads include dir under Windows
#-----------------------------------------------------------------------------

if (WIN32)
    if (NOT COIN_THREADS_INC_PATH STREQUAL "None")
        include_directories(${COIN_THREADS_INC_PATH})
    endif ()
    if (NOT COIN_THREADS_LIB_PATH STREQUAL "None")
        link_directories(${COIN_THREADS_LIB_PATH})
    endif ()
endif ()

# Check for a fortran compiler
include(CMakeDetermineFortranCompiler)
if (NOT CMAKE_Fortran_COMPILER)
    message(STATUS "WARNING: fortran compiler not found. Disabling f77/f95 bindings")
endif ()

# Define FORTRAN_INTEGER_TYPE for Ipopt.
set(FORTRAN_INTEGER_TYPE int)

#-----------------------------------------------------------------------------
# Detect name mangling convention used between Fortran and C
#-----------------------------------------------------------------------------

# TODO: this fails for flang with gcc & g++

if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran)

    safe_tmp(_LOCAL_TMP_F_FLAGS CMAKE_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-lto")

    include(FortranCInterface)

    FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/F77Mangle.h
            MACRO_NAMESPACE "F77_"
            SYMBOL_NAMESPACE "F77_")

    file(STRINGS ${CMAKE_CURRENT_BINARY_DIR}/F77Mangle.h CONTENTS REGEX "F77_GLOBAL\\(.*,.*\\) +(.*)")
    string(REGEX MATCH "F77_GLOBAL\\(.*,.*\\) +(.*)" RESULT ${CONTENTS})
    set(F77_FUNC         "F77_FUNC(name,NAME) ${CMAKE_MATCH_1}")
    set(COIN_LAPACK_FUNC "COIN_LAPACK_FUNC(name,NAME) ${CMAKE_MATCH_1}")

    file(STRINGS ${CMAKE_CURRENT_BINARY_DIR}/F77Mangle.h CONTENTS REGEX "F77_GLOBAL_\\(.*,.*\\) +(.*)")
    string(REGEX MATCH "F77_GLOBAL_\\(.*,.*\\) +(.*)" RESULT ${CONTENTS})
    set(F77_FUNC_         "F77_FUNC_(name,NAME) ${CMAKE_MATCH_1}")
    set(COIN_LAPACK_FUNC_ "COIN_LAPACK_FUNC_(name,NAME) ${CMAKE_MATCH_1}")
    restore_tmp(CMAKE_Fortran_FLAGS _LOCAL_TMP_F_FLAGS)
else ()
    set(F77_FUNC          "F77_FUNC(name,NAME)  name##_")
    set(F77_FUNC_         "F77_FUNC_(name,NAME) name##__")
    set(COIN_LAPACK_FUNC  "COIN_LAPACK_FUNC(name,NAME)  name##_")
    set(COIN_LAPACK_FUNC_ "COIN_LAPACK_FUNC_(name,NAME) name##__")
endif ()


set(F77_DUMMY_MAIN "" CACHE STRING "Define to dummy 'main' function (if any) required to link to the Fortran libraries.")
set(FC_DUMMY_MAIN  "" CACHE STRING "Define to dummy 'main' function (if any) required to link to the Fortran libraries.")
option(FC_DUMMY_MAIN_EQ_F77 "Define if F77 and FC dummy 'main' functions are identical." OFF)

if (FC_DUMMY_MAIN_EQ_F77)
    set(FC_DUMMY_MAIN "${F77_DUMMY_MAIN}")
endif ()

mark_as_advanced(F77_FUNC
        F77_FUNC_
        F77_DUMMY_MAIN
        FC_DUMMY_MAIN
        FC_DUMMY_MAIN_EQ_F77)

# End of coverage

if (USE_FAST_CODE)
    set(CMAKE_CXX_FLAGS       "${CMAKE_CXX_FLAGS} -O3 -fomit-frame-pointer ")
    set(CMAKE_C_FLAGS         "${CMAKE_C_FLAGS} -O3 -fomit-frame-pointer ")
    set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -O3 ")

    if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno -fp-trap=none ")
        set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fno-math-errno -fp-trap=none ")
    endif ()

    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -fno-math-errno -fno-trapping-math ")
        set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -ffast-math -fno-math-errno -fno-trapping-math ")
    endif ()

endif ()

add_definitions(-DHAVE_CONFIG_H)

add_subdirectory(ThirdParty)

#
# HSL Management
#

# TODO: add if we support HSL from source build
#
# if (HSL_HAS_METIS)
#     add_definitions(-DCOINHSL_HAS_METIS)
# endif ()
#

add_subdirectory(Ipopt)

#
# Packaging
#

set(CPACK_PACKAGE_NAME      "${PACKAGE_NAME}")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${IPOPT_VERSION}-${CMAKE_SYSTEM_NAME}")

message(STATUS "Package filename: ${CPACK_PACKAGE_FILE_NAME}")

set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PACKAGE_NAME}")
set(CPACK_PACKAGE_VERSION "${IPOPT_VERSION}")

set(CPACK_PACKAGE_RELOCATABLE TRUE)

if (WIN32)
    set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME}")
    set(CPACK_MONOLITHIC_INSTALL ON)
    set(CPACK_NSIS_CONTACT "${PACKAGE_BUGREPORT}")
    set(CPACK_NSIS_MODIFY_PATH ON)
    set(CPACK_NSIS_PACKAGE_NAME "Ipopt ${IPOPT_VERSION}")
elseif (APPLE)
    set(CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE ON)
    set(CPACK_PACKAGE_DEFAULT_LOCATION "/opt/${CPACK_PACKAGE_NAME}")
    set(CPACK_PACKAGING_INSTALL_PREFIX "/")
else (WIN32)
    set(CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE ON)
    set(CPACK_PACKAGE_DEFAULT_LOCATION "/opt/${CPACK_PACKAGE_NAME}")
    set(CPACK_PACKAGING_INSTALL_PREFIX "/")

    set(CPACK_DEBIAN_PACKAGE_SECTION "Libraries")
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "${PACKAGE_BUGREPORT}")

    set(CPACK_RPM_COMPONENT_INSTALL ON)
    set(CPACK_RPM_PACKAGE_RELOCATABLE ON)
    set(CPACK_RPM_PACKAGE_LICENSE "Copyrighted")
    set(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
    set(CPACK_RPM_PACKAGE_URL "${PACKAGE_URL}")
    set(CPACK_RPM_PACKAGE_SUMMARY "Ipopt is a non linear constrained solver library from CoinOR.")
    set(CPACK_RPM_PACKAGE_DESCRIPTION
            "The CoinUtils suite includes:

    * The Ipopt non linear constrained solver.
"
            )
endif ()

include(CPack)
include(CTest)

message(STATUS "Installation directory: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

get_directory_property(ALL_INCLUDES INCLUDE_DIRECTORIES)
get_directory_property(ALL_LINKS    LINK_DIRECTORIES)
get_directory_property(ALL_DEFS     COMPILE_DEFINITIONS)

message(STATUS "Include directories:${color_reset}     ${ALL_INCLUDES}")
message(STATUS "Link directories:${color_reset}        ${ALL_LINKS}")
message(STATUS "Compilation definitions:${color_reset} ${ALL_DEFS}")

#
# Install part
#

# install(FILES ${Ipopt_DIR}/LICENSE ${Ipopt_DIR}/README.md
#         DESTINATION share/coin/doc/Ipopt/)

# TODO: remove this temporary solution
restore_tmp(CMAKE_C_FLAGS        _GLOBL_TMP_CFLAGS       )
restore_tmp(CMAKE_CXX_FLAGS      _GLOBL_TMP_CXXFLAGS     )
restore_tmp(CMAKE_Fortran_FLAGS  _GLOBL_TMP_FortranFLAGS )
