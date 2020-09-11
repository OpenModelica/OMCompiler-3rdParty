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
set(enable_java_finalization OFF CACHE BOOL "Support for java finalization")
set(enable_gcj_support OFF CACHE BOOL "Support for gcj")
set(enable_large_config ON CACHE BOOL "Optimize for large heap or root set")
omc_add_subdirectory(gc)
add_library(omc::3rd::gc ALIAS gc-lib)
# make sure every target that links to gc-lib has its sources
# compiled with -DGC_THREADS 
target_compile_definitions(gc-lib INTERFACE GC_THREADS)

# On Windows with MinGW OM uses pthreads-win32.
# ideally this should made in gc's own cmakelists.txt.
# but that file is rather messy. Until we either get a proper cmake
# implementation from upstream or we make one ourselves, we do it
# here. I also don't want to make it require recent cmake just for this line.
if(MINGW)
    target_compile_definitions(gc-lib PUBLIC GC_WIN32_PTHREADS)
endif(MINGW)




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
omc_add_subdirectory(sundials)


