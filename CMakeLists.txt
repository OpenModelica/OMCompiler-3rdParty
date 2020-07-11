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
target_include_directories(fmilib_shared INTERFACE ${FMILibrary_SOURCE_DIR}/install/include)
add_library(omc::3rd::fmilib::shared ALIAS fmilib_shared)
add_library(omc::3rd::FMIL::minizip ALIAS minizip)



# Boehm GC
set(CMAKE_USE_WIN32_THREADS_INIT ON)
set(enable_java_finalization OFF)
set(enable_gcj_support OFF)
set(enable_large_config ON)
omc_add_subdirectory(gc)
add_library(omc::3rd::gc ALIAS gc-lib)



# graphstream
omc_add_subdirectory(graphstream)
add_library(omc::3rd::netstream ALIAS netstream)



# libzmq
set (ZMQ_BUILD_TESTS OFF CACHE BOOL "Build the tests for ZeroMQ")
omc_add_subdirectory(libzmq)
add_library(omc::3rd::libzmq ALIAS libzmq)



# lpsolve55
omc_add_subdirectory(lpsolve)
add_library(omc::3rd::lpsolve55 ALIAS lpsolve55)



# metis
omc_add_subdirectory(metis-5.1.0)
add_library(omc::3rd::metis ALIAS metis)



# Intel TBB
set(TBB_BUILD_TESTS OFF CACHE BOOL "Build TBB tests and enable testing infrastructure")
set(TBB_BUILD_SHARED OFF CACHE BOOL "Build TBB shared library")
set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL "Build TBB malloc library")
set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL "Build TBB malloc proxy library")
omc_add_subdirectory(tbb_cmake)
add_library(omc::3rd::tbb::static ALIAS tbb_static)
# add_library(omc::3rd::tbb::shared ALIAS tbb)



# # regex
# omc_add_subdirectory(regex-0.12)

# SuiteSparse
omc_add_subdirectory(SuiteSparse)

# sundials
omc_add_subdirectory(sundials)


