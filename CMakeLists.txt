project(OMCompiler_3rdParty)



# antlr
omc_add_subdirectory(antlr)
add_library(omc::3rd::omantlr3 ALIAS omantlr3)



# CMinpack
omc_add_subdirectory(CMinpack)



# # cppzmq
# omc_add_subdirectory(cppzmq)



# FMIL
set (FMILIB_GENERATE_DOXYGEN_DOC OFF CACHE BOOL "Generate doxygen doc target")
set (FMILIB_BUILD_TESTS OFF CACHE BOOL "Build tests")
omc_add_subdirectory(FMIL)
omc_add_to_report(FMILIB_FMI_PLATFORM)
add_library(omc::3rd::fmilib::shared ALIAS fmilib_shared)
add_library(omc::3rd::FMIL::minizip ALIAS minizip)



# Boehm GC
set(CMAKE_USE_WIN32_THREADS_INIT ON)
set(enable_java_finalization OFF)
set(enable_gcj_support OFF)
set(enable_large_config ON)
omc_add_subdirectory(gc)
add_library(omc::3rd::gc ALIAS gc-lib)



# libzmq
set (ZMQ_BUILD_TESTS OFF CACHE BOOL "Build the tests for ZeroMQ")
omc_add_subdirectory(libzmq)
add_library(omc::3rd::libzmq ALIAS libzmq)



# lpsolve55
omc_add_subdirectory(lpsolve)
add_library(omc::3rd::lpsolve55 ALIAS lpsolve55)



# metis
omc_add_subdirectory(metis-5.1.0)



# Intel TBB
set(TBB_BUILD_TESTS OFF CACHE BOOL "Build TBB tests and enable testing infrastructure")
set(TBB_BUILD_SHARED OFF CACHE BOOL "Build TBB shared library")
set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL "Build TBB malloc library")
set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL "Build TBB malloc proxy library")
omc_add_subdirectory(tbb_cmake)



# # regex
# omc_add_subdirectory(regex-0.12)

# SuiteSparse
omc_add_subdirectory(SuiteSparse)

# sundials
omc_add_subdirectory(sundials)


