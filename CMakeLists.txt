project(OMCompiler_3rdParty)



# antlr
omc_add_subdirectory(antlr)
add_library(omc::3rd::omantlr3 ALIAS omantlr3)



# CMinpack
omc_add_subdirectory(CMinpack)
add_library(omc::3rd::cminpack ALIAS cminpack)



# # cppzmq
# omc_add_subdirectory(cppzmq)



# FMIL
set (FMILIB_GENERATE_DOXYGEN_DOC OFF CACHE BOOL "Generate doxygen doc target")
set (FMILIB_BUILD_TESTS OFF CACHE BOOL "Build tests")
omc_add_subdirectory(FMIL)
# For now we deal with FMIL's nonsensical structure here. They really need to
# fix up their cmake usage. It is their default build system and it is used in
# very contrived and unusual ways.
target_include_directories(fmilib INTERFACE ${FMILibrary_SOURCE_DIR}/install/include)
add_library(omc::3rd::fmilib::static ALIAS fmilib)
add_library(omc::3rd::FMIL::minizip ALIAS minizip)



# Boehm GC
# We use pthreads API even on Windows
set(CMAKE_USE_PTHREADS_INIT ON)

set(GC_BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")
set(enable_java_finalization OFF CACHE BOOL "Support for java finalization")
set(enable_gcj_support OFF CACHE BOOL "Support for gcj")
set(enable_large_config ON CACHE BOOL "Optimize for large heap or root set")
omc_add_subdirectory(gc)
target_include_directories(omcgc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/gc/include)
# make sure every target that links to gc-lib has its sources
# compiled with -DGC_WIN32_PTHREADS (for pthreads on Windows, i.e., our MingW)
# Or -DGC_THREADS (for auto detection on other systems.)
# On Windows with MinGW OM uses pthreads-win32. GC_WIN32_PTHREADS is required
# to be set explicitly for use of pthreads API on Windows.
if(MINGW)
    target_compile_definitions(omcgc PUBLIC GC_WIN32_PTHREADS)
else(MINGW)
    target_compile_definitions(omcgc PUBLIC GC_THREADS)
endif(MINGW)

# Finally add an alias for clarity purposes.
add_library(omc::3rd::omcgc ALIAS omcgc)





# graphstream
omc_add_subdirectory(graphstream)
add_library(omc::3rd::netstream ALIAS netstream)



# libzmq
set (ZMQ_BUILD_TESTS OFF CACHE BOOL "Build the tests for ZeroMQ")
set (ZMQ_BUILD_STATIC ON CACHE BOOL "Whether or not to build the static object")
set (ZMQ_BUILD_SHARED OFF CACHE BOOL "Whether or not to build the shared object")
omc_add_subdirectory(libzmq)
add_library(omc::3rd::libzmq ALIAS libzmq-static)



# lpsolve55
omc_add_subdirectory(lpsolve)
add_library(omc::3rd::lpsolve55 ALIAS lpsolve55)



# metis
omc_add_subdirectory(metis-5.1.0)
add_library(omc::3rd::metis ALIAS metis)
target_include_directories(metis INTERFACE metis-5.1.0/include)


# ModelicaExternalC
omc_add_subdirectory(ModelicaExternalC)
add_library(omc::3rd::modelica_external_c ALIAS ModelicaExternalC)
add_library(omc::3rd::modelica_mat_io ALIAS ModelicaMatIO)
add_library(omc::3rd::modelica_io ALIAS ModelicaIO)


# Intel TBB
# set(TBB_BUILD_TESTS OFF CACHE BOOL "Build TBB tests and enable testing infrastructure")
# set(TBB_BUILD_SHARED OFF CACHE BOOL "Build TBB shared library")
# set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL "Build TBB malloc library")
# set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL "Build TBB malloc proxy library")
# omc_add_subdirectory(tbb_cmake)
# add_library(omc::3rd::tbb::static ALIAS tbb_static)
# add_library(omc::3rd::tbb::shared ALIAS tbb)



# # regex
# omc_add_subdirectory(regex-0.12)

# SuiteSparse
omc_add_subdirectory(SuiteSparse)

# sundials
omc_add_subdirectory(sundials-5.4.0)


