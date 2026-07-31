project(OMCompiler_3rdParty)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ryu
omc_add_subdirectory(ryu)
add_library(omc::3rd::ryu ALIAS ryu)

# zlib
option(OM_USE_SYSTEM_ZLIB "Use system zlib" OFF)
if (OM_USE_SYSTEM_ZLIB)
  set(ZLIB_USE_STATIC_LIBS ON)
  find_package(ZLIB MODULE REQUIRED)
  unset(ZLIB_USE_STATIC_LIBS)
  add_library(omc::3rd::zlib INTERFACE)
  target_link_libraries(omc::3rd::zlib INTERFACE ZLIB::ZLIB)
else ()
  omc_add_subdirectory(zlib)
  add_library(omc::3rd::zlib ALIAS zlib)
endif ()

# antlr
omc_add_subdirectory(antlr)
add_library(omc::3rd::omantlr3 ALIAS omantlr3)

# CDaskr
omc_add_subdirectory(Cdaskr)
add_library(omc::3rd::cdaskr ALIAS cdaskr)

# CMinpack
option(OM_USE_SYSTEM_CMINPACK "Use system cminpack" OFF)
if (OM_USE_SYSTEM_CMINPACK)
  find_package(CMinpack CONFIG REQUIRED)
  add_library(cminpack INTERFACE)
  target_link_libraries(cminpack INTERFACE cminpack::cminpack)
else ()
  omc_add_subdirectory(CMinpack)
endif ()
add_library(omc::3rd::cminpack ALIAS cminpack)

# # cppzmq
# omc_add_subdirectory(cppzmq)



# FMIL
set (FMILIB_GENERATE_DOXYGEN_DOC OFF CACHE BOOL "Generate doxygen doc target")
set (FMILIB_BUILD_TESTS OFF CACHE BOOL "Build tests")
set (FMILIB_BUILD_SHARED_LIB OFF CACHE BOOL "Build the library as shared (dll/so/dylib).")
# Build the bundled Expat as a subdirectory. FMIL builds it as an ExternalProject by
# default, which leaves an imported target that the rest of the build cannot link to.
set (FMILIB_EXPAT_AS_SUBPROJECT ON CACHE BOOL "Build the bundled Expat with add_subdirectory")
# zlib is added above, so FMIL picks up that target instead of building its own copy.
omc_add_subdirectory(FMIL)

# For now we deal with FMIL's nonsensical structure here. They really need to
# fix up their structure as a whole! It is quite non-conventional.
# They also need to fix their cmake usage. It is their default build system and it is used in
# very contrived and unusual ways.

# We create a top level 'include' directory that matches their include structure when FMIL is installed-to-be-used.
# This is how every library should be. That way when you install the library you just change the include
# path and every include will be resolved as it was during build time.  Ideally it should even be in  'include/FMIL/'
# but that won't match how they install their files right now.
file(MAKE_DIRECTORY ${FMILibrary_SOURCE_DIR}/include)
# Do not ask me why the fmilib.h is in the cmake.config directory.
file(COPY ${FMILibrary_SOURCE_DIR}/Config.cmake/fmilib.h DESTINATION ${FMILibrary_SOURCE_DIR}/include)

file(COPY ${FMILibrary_SOURCE_DIR}/src/Import/include/FMI DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/src/Import/include/FMI1 DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/src/Import/include/FMI2 DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/src/Import/include/FMI3 DESTINATION ${FMILibrary_SOURCE_DIR}/include)

file(COPY ${FMILibrary_SOURCE_DIR}/src/Util/include/FMI DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/src/Util/include/FMI1 DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/src/Util/include/FMI2 DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/src/Util/include/FMI3 DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/src/Util/include/JM DESTINATION ${FMILibrary_SOURCE_DIR}/include)

file(COPY ${FMILibrary_SOURCE_DIR}/ThirdParty/FMI/default/FMI1 DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/ThirdParty/FMI/default/FMI2 DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/ThirdParty/FMI/default/FMI3 DESTINATION ${FMILibrary_SOURCE_DIR}/include)

file(COPY ${FMILibrary_BINARY_DIR}/fmilib_config_include/fmilib_config.h DESTINATION ${FMILibrary_SOURCE_DIR}/include)

# We give this new directory as include dir for targets that depend on FMIL.
target_include_directories(fmilib INTERFACE ${FMILibrary_SOURCE_DIR}/include)
add_library(omc::3rd::fmilib ALIAS fmilib)
add_library(omc::3rd::FMIL::minizip ALIAS minizip)
add_library(omc::3rd::FMIL::expat ALIAS expat)



# Boehm GC
# We use pthreads API even on Windows
set(CMAKE_USE_PTHREADS_INIT ON)

option(GC_BUILD_SHARED_LIBS "Build shared libraries" ON)
option(enable_java_finalization "Support for java finalization" OFF)
option(enable_gcj_support "Support for gcj" OFF)
option(enable_large_config "Optimize for large heap or root set" ON)
omc_add_subdirectory(gc)
target_include_directories(omcgc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/gc/include)
# make sure every target that links to gc-lib has its sources
# compiled with -DGC_WIN32_PTHREADS (for pthreads on Windows, i.e., OMDev or MSVC)
# Or -DGC_THREADS (for auto detection on other systems.)
# Even on Windows OM uses pthreads. GC_WIN32_PTHREADS is required
# to be set explicitly for use of pthreads API on Windows.
if(WIN32)
    target_compile_definitions(omcgc PUBLIC GC_WIN32_PTHREADS)
else()
    target_compile_definitions(omcgc PUBLIC GC_THREADS)
endif()

# Finally add an alias for clarity purposes.
add_library(omc::3rd::omcgc ALIAS omcgc)





# libffi
option(OM_USE_SYSTEM_LIBFFI "Use system libffi" OFF)
if (OM_USE_SYSTEM_LIBFFI)
  find_package(LibFFI MODULE REQUIRED)
  add_library(ffi INTERFACE)
  target_include_directories(ffi INTERFACE ${LIBFFI_INCLUDE_DIRS})
  target_link_libraries(ffi INTERFACE ${LIBFFI_LIBRARIES})
elseif (CMAKE_CROSSCOMPILING)
  add_library(ffi INTERFACE)
else ()
  omc_add_subdirectory(libffi)
  target_include_directories(ffi INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/libffi/include)
  # Currently the libffi setup puts the configure-generated headers in the
  # `build directory`/include folder
  target_include_directories(ffi INTERFACE ${libffi_BINARY_DIR}/include)
endif ()
add_library(omc::3rd::ffi ALIAS ffi)




# libzmq
option(OM_USE_SYSTEM_ZMQ "Use system zmq" OFF)
if (OM_USE_SYSTEM_ZMQ)
  find_package(ZeroMQ CONFIG REQUIRED)
  add_library(zmq INTERFACE)
  target_link_libraries(zmq INTERFACE libzmq)
  add_library(omc::3rd::libzmq ALIAS zmq)
else ()
  set (ZMQ_BUILD_TESTS OFF CACHE BOOL "Build the tests for ZeroMQ")
  set (ZMQ_BUILD_STATIC ON CACHE BOOL "Whether or not to build the static object")
  set (ZMQ_BUILD_SHARED OFF CACHE BOOL "Whether or not to build the shared object")
  omc_add_subdirectory(libzmq)
  add_library(omc::3rd::libzmq ALIAS libzmq-static)
endif ()


# simdjson
omc_add_subdirectory(simdjson)
add_library(omc::3rd::simdjson ALIAS simdjson)

#lis
omc_add_subdirectory(lis-1.4.12)
add_library(omc::3rd::lis ALIAS lis)


# metis
omc_add_subdirectory(metis-5.1.0)
add_library(omc::3rd::metis ALIAS metis)
target_include_directories(metis INTERFACE metis-5.1.0/include)

omc_add_subdirectory(open62541)
add_library(omc::3rd::opcua ALIAS opcua)

# Intel oneTBB (github.com/uxlfoundation/oneTBB)
set(TBB_TEST OFF CACHE BOOL "Enable oneTBB testing" FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL "Enable oneTBB examples" FORCE)
set(TBB_STRICT OFF CACHE BOOL "Treat oneTBB compiler warnings as errors" FORCE)
set(TBBMALLOC_BUILD OFF CACHE BOOL "Enable tbbmalloc build" FORCE)
set(TBB_INSTALL OFF CACHE BOOL "Enable oneTBB installation" FORCE)
# Build oneTBB statically, matching the previous (wjakob/tbb 2020.2) setup.
# Scope BUILD_SHARED_LIBS to the tbb subdirectory only so other 3rdParty libs are unaffected.
set(_omc_saved_build_shared_libs ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF)
omc_add_subdirectory(oneTBB)
set(BUILD_SHARED_LIBS ${_omc_saved_build_shared_libs})
add_library(omc::3rd::tbb ALIAS tbb)

# Install oneTBB (under its own default library name) next to the other omc
# runtime libraries so generated --parmodauto simulators can link against it.
install(TARGETS tbb ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Resolve oneTBB's default library name for the active build type (it appends a
# _debug postfix in Debug builds) and expose it so the code generator can emit the
# matching -l flag via Autoconf.parModelicaAutoLibs. Unix Makefiles/Ninja are
# single-config, so CMAKE_BUILD_TYPE selects the right postfix here.
get_target_property(_omc_tbb_output_name tbb OUTPUT_NAME)
if(NOT _omc_tbb_output_name)
  set(_omc_tbb_output_name tbb)
endif()
set(_omc_tbb_postfix "")
if(CMAKE_BUILD_TYPE)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _omc_tbb_cfg)
  get_target_property(_omc_tbb_cfg_postfix tbb ${_omc_tbb_cfg}_POSTFIX)
  if(_omc_tbb_cfg_postfix)
    set(_omc_tbb_postfix "${_omc_tbb_cfg_postfix}")
  endif()
endif()
set(OMC_PARMODELICA_TBB_LINK_NAME "${_omc_tbb_output_name}${_omc_tbb_postfix}"
    CACHE INTERNAL "oneTBB library name generated --parmodauto simulators link against")


# regex
if(MSVC)
  omc_add_subdirectory(regex-0.12)
  add_library(omc::3rd::regex ALIAS omcregex)
endif()


# SuiteSparse v7 (DrTimothyAldenDavis/SuiteSparse, submodule pinned to v7.12.2)
# The hand-written 5.8.1 wrapper was replaced with a submodule and we delegate
# the build to upstream by including the top-level SuiteSparse/CMakeLists.txt.
# Only the subset required by OMC and Sundials is enabled.
set(SUITESPARSE_ENABLE_PROJECTS "suitesparse_config;amd;colamd;btf;klu;umfpack" CACHE STRING "SuiteSparse projects to build" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_STATIC_LIBS ON CACHE BOOL "" FORCE)
set(SUITESPARSE_USE_OPENMP OFF CACHE BOOL "" FORCE)
set(SUITESPARSE_USE_CUDA OFF CACHE BOOL "" FORCE)
set(SUITESPARSE_USE_FORTRAN OFF CACHE BOOL "" FORCE)
set(SUITESPARSE_DEMOS OFF CACHE BOOL "" FORCE)
set(SUITESPARSE_USE_STRICT OFF CACHE BOOL "" FORCE)
set(KLU_USE_CHOLMOD OFF CACHE BOOL "" FORCE)
set(UMFPACK_USE_CHOLMOD OFF CACHE BOOL "" FORCE)

omc_add_subdirectory(SuiteSparse)

add_library(omc::3rd::suitesparse::config ALIAS SuiteSparseConfig_static)
add_library(omc::3rd::suitesparse::amd ALIAS AMD_static)
add_library(omc::3rd::suitesparse::colamd ALIAS COLAMD_static)
add_library(omc::3rd::suitesparse::btf ALIAS BTF_static)
add_library(omc::3rd::suitesparse::klu ALIAS KLU_static)
add_library(omc::3rd::suitesparse::umfpack ALIAS UMFPACK_static)


# SUNDIALS (LLNL/sundials, submodule pinned to v7.8.0)
# The patched 5.4.0 source copy was replaced by a pristine submodule, so
# everything we used to patch into those sources is configured from here.
# See sundials.txt for what those patches were and why they are (or are no
# longer) needed.

# Static only, like the rest of 3rdParty. BUILD_STATIC_LIBS/BUILD_SHARED_LIBS are
# already forced above for SuiteSparse; SUNDIALS reads the same two variables.

# The OMC runtimes and the generated simulation code assume doubles and a 64 bit
# sunindextype (which also has to match SuiteSparse's int64 KLU entry points).
set(SUNDIALS_PRECISION "DOUBLE" CACHE STRING "SUNDIALS floating point precision" FORCE)
set(SUNDIALS_INDEX_SIZE "64" CACHE STRING "SUNDIALS sunindextype size in bits" FORCE)

# KLU. SUNDIALS' FindKLU prefers an already existing SuiteSparse::KLU target as
# long as we do not point KLU_INCLUDE_DIR/KLU_LIBRARY_DIR/KLU_LIBRARY somewhere
# else, so it picks up the in-tree SuiteSparse built above on its own - do not
# hand it anything. The compatibility checks have to be off: they try_compile()
# against SUNDIALS::KLU, which here is an alias for an in-tree target rather than
# an IMPORTED one, so it is not visible to the separate try_compile project.
# NOTE: FindKLU tries a find_package(KLU CONFIG) first. On distros that ship
# SuiteSparse' CMake package files this prints a "found KLUConfig.cmake but it set
# KLU_FOUND to FALSE" warning; it is harmless. Because SuiteSparse::KLU is already
# a target by the time we get here, the system package can neither be imported nor
# override it, so the in-tree KLU_static wins either way.
set(SUNDIALS_ENABLE_KLU ON CACHE BOOL "Enable KLU support" FORCE)
set(SUNDIALS_ENABLE_KLU_CHECKS OFF CACHE BOOL "Enable KLU compatibility checks" FORCE)

set(SUNDIALS_ENABLE_LAPACK ON CACHE BOOL "Enable LAPACK support" FORCE)
# SUNDIALS needs a Fortran compiler to infer the LAPACK name-mangling scheme. If
# OpenModelica is built without Fortran, spell the scheme out instead - almost
# every Fortran-to-C converted LAPACK we have seen is lower case with one
# underscore.
if(NOT OM_OMC_ENABLE_FORTRAN)
  set(SUNDIALS_LAPACK_CASE "lower" CACHE STRING "Case of LAPACK function names" FORCE)
  set(SUNDIALS_LAPACK_UNDERSCORES "one" CACHE STRING "Underscores appended to LAPACK function names" FORCE)
endif()

# We use the C interfaces only, and want neither examples, benchmarks nor tests.
set(SUNDIALS_ENABLE_FORTRAN OFF CACHE BOOL "Enable SUNDIALS Fortran interfaces" FORCE)
set(SUNDIALS_ENABLE_C_EXAMPLES OFF CACHE BOOL "Build SUNDIALS C examples" FORCE)
set(SUNDIALS_ENABLE_CXX_EXAMPLES OFF CACHE BOOL "Build SUNDIALS C++ examples" FORCE)
set(SUNDIALS_ENABLE_EXAMPLES_INSTALL OFF CACHE BOOL "Install SUNDIALS examples" FORCE)
set(SUNDIALS_ENABLE_BENCHMARKS OFF CACHE BOOL "Enable the SUNDIALS benchmark suite" FORCE)
set(SUNDIALS_TEST_ENABLE_UNIT_TESTS OFF CACHE BOOL "Enable SUNDIALS unit tests" FORCE)

# Upstream turns the extensive error checks on by default for Debug builds (and
# warns about the performance hit while configuring). Keep SUNDIALS behaving the
# same way in all OMC build types.
set(SUNDIALS_ENABLE_ERROR_CHECKS OFF CACHE BOOL "Enable SUNDIALS error checking" FORCE)

# OMC installs the SUNDIALS headers under <includedir>/sundials, which
# SimulationRuntime/cpp, SimulationRuntime/OMSICpp and
# SimulationRuntime/c/Makefile.common all hardcode. SUNDIALS installs everything
# relative to CMAKE_INSTALL_INCLUDEDIR, so bump that for the subdirectory only.
set(_omc_saved_install_includedir "${CMAKE_INSTALL_INCLUDEDIR}")
set(CMAKE_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/sundials")

omc_add_subdirectory(sundials)

set(CMAKE_INSTALL_INCLUDEDIR "${_omc_saved_install_includedir}")
unset(_omc_saved_install_includedir)

# SUNDIALS declares CMAKE_C_STANDARD and CMAKE_C_EXTENSIONS as *cache* variables
# (cmake/SundialsSetupCompilers.cmake), so configuring it would otherwise drag the
# whole OpenModelica build down to C99. OMC deliberately leaves the C standard at
# the compiler default, see the top-level CMakeLists.txt. SUNDIALS' own targets
# already captured the standard when they were created above, so dropping the
# cache entries here only affects everything that comes after.
unset(CMAKE_C_STANDARD CACHE)
unset(CMAKE_C_EXTENSIONS CACHE)

## Add aliases for the static libs. For readability.
## Linking one of these already provides the SUNDIALS include directories and the
## SUNDIALS_STATIC_DEFINE compile definition transitively - sundials_add_library()
## attaches both as PUBLIC usage requirements - so unlike with 5.4.0 we do not
## need an interface library of our own on top.
add_library(omc::3rd::sundials::cvode ALIAS sundials_cvode_static)
add_library(omc::3rd::sundials::ida ALIAS sundials_ida_static)
add_library(omc::3rd::sundials::idas ALIAS sundials_idas_static)
add_library(omc::3rd::sundials::kinsol ALIAS sundials_kinsol_static)
add_library(omc::3rd::sundials::sunlinsolklu ALIAS sundials_sunlinsolklu_static)
add_library(omc::3rd::sundials::sunlinsollapackdense ALIAS sundials_sunlinsollapackdense_static)

if(OM_OMC_ENABLE_OPTIMIZATION)
    add_library(EXTERNAL_METIS_TARGET ALIAS metis)
    omc_add_subdirectory(moo-0.1.0)
    add_library(omc::3rd::moo ALIAS moo)
    add_library(omc::3rd::ipopt ALIAS ipopt)
endif() # OM_OMC_ENABLE_OPTIMIZATION

if(OM_OMC_ENABLE_PRIMME)
    omc_add_subdirectory(primme-3.2.3)
    add_library(omc::3rd::primme ALIAS primme)
endif() # OM_OMC_ENABLE_PRIMME

if(OM_OMC_ENABLE_COLPACK)
    omc_add_subdirectory(ColPack)
    add_library(omc::3rd::colpack ALIAS colpack)
endif() # OM_OMC_ENABLE_COLPACK
