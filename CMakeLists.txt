project(OMCompiler_3rdParty)

# antlr
omc_add_subdirectory(zlib)
add_library(omc::3rd::zlib ALIAS zlib)

# antlr
omc_add_subdirectory(antlr)
add_library(omc::3rd::omantlr3 ALIAS omantlr3)

# CDaskr
omc_add_subdirectory(Cdaskr)
add_library(omc::3rd::cdaskr ALIAS cdaskr)

# CMinpack
omc_add_subdirectory(CMinpack)
add_library(omc::3rd::cminpack ALIAS cminpack)


# # cppzmq
# omc_add_subdirectory(cppzmq)



# FMIL
set (FMILIB_GENERATE_DOXYGEN_DOC OFF CACHE BOOL "Generate doxygen doc target")
set (FMILIB_BUILD_TESTS OFF CACHE BOOL "Build tests")
set (FMILIB_BUILD_SHARED_LIB OFF CACHE BOOL "Build the library as shared (dll/so/dylib).")
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

file(COPY ${FMILibrary_SOURCE_DIR}/src/Util/include/FMI DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/src/Util/include/FMI1 DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/src/Util/include/FMI2 DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/src/Util/include/JM DESTINATION ${FMILibrary_SOURCE_DIR}/include)

file(COPY ${FMILibrary_SOURCE_DIR}/ThirdParty/FMI/default/FMI1 DESTINATION ${FMILibrary_SOURCE_DIR}/include)
file(COPY ${FMILibrary_SOURCE_DIR}/ThirdParty/FMI/default/FMI2 DESTINATION ${FMILibrary_SOURCE_DIR}/include)

file(COPY ${FMILibrary_BINARY_DIR}/fmilib_config.h DESTINATION ${FMILibrary_SOURCE_DIR}/include)

# We give this new directory as include dir for targets that depend on FMIL.
target_include_directories(fmilib INTERFACE ${FMILibrary_SOURCE_DIR}/include)
add_library(omc::3rd::fmilib ALIAS fmilib)
add_library(omc::3rd::FMIL::minizip ALIAS minizip)
add_library(omc::3rd::FMIL::expat ALIAS expat)



# Boehm GC
# We use pthreads API even on Windows
set(CMAKE_USE_PTHREADS_INIT ON)

option(GC_BUILD_SHARED_LIBS "Build shared libraries" OFF)
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





#libffi
omc_add_subdirectory(libffi)
target_include_directories(ffi INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/libffi_cmake/include)
# Currently the libffi setup puts the configure-generated headers in the
# `build directory`/include folder
target_include_directories(ffi INTERFACE ${libffi_BINARY_DIR}/include)
add_library(omc::3rd::ffi ALIAS ffi)




# libzmq
set (ZMQ_BUILD_TESTS OFF CACHE BOOL "Build the tests for ZeroMQ")
set (ZMQ_BUILD_STATIC ON CACHE BOOL "Whether or not to build the static object")
set (ZMQ_BUILD_SHARED OFF CACHE BOOL "Whether or not to build the shared object")
omc_add_subdirectory(libzmq)
add_library(omc::3rd::libzmq ALIAS libzmq-static)



# lpsolve55
if(OM_OMC_USE_LPSOLVE AND OM_OMC_BUILD_LPSOLVE)
  omc_add_subdirectory(lpsolve)
  add_library(omc::3rd::lpsolve55 ALIAS lpsolve55)
endif()

#lis
omc_add_subdirectory(lis-1.4.12)
add_library(omc::3rd::lis ALIAS lis)


# metis
omc_add_subdirectory(metis-5.1.0)
add_library(omc::3rd::metis ALIAS metis)
target_include_directories(metis INTERFACE metis-5.1.0/include)


omc_add_subdirectory(open62541)
add_library(omc::3rd::opcua ALIAS opcua)

# Intel TBB
option(TBB_BUILD_TESTS "Build TBB tests and enable testing infrastructure" OFF)
option(TBB_BUILD_SHARED "Build TBB shared library" OFF)
option(TBB_BUILD_TBBMALLOC "Build TBB malloc library" OFF)
option(TBB_BUILD_TBBMALLOC_PROXY "Build TBB malloc proxy library" OFF)
omc_add_subdirectory(tbb)
add_library(omc::3rd::tbb ALIAS tbb_static)
# add_library(omc::3rd::tbb::shared ALIAS tbb)



# # regex
# omc_add_subdirectory(regex-0.12)

# SuiteSparse
omc_add_subdirectory(SuiteSparse-5.8.1)
add_library(omc::3rd::suitesparse::klu ALIAS klu)
target_include_directories(klu INTERFACE SuiteSparse-5.8.1/KLU/Include)

add_library(omc::3rd::suitesparse::umfpack ALIAS umfpack)
target_include_directories(umfpack INTERFACE SuiteSparse-5.8.1/UMFPACK/Include)

add_library(omc::3rd::suitesparse::amd ALIAS amd)
target_include_directories(amd INTERFACE SuiteSparse-5.8.1/AMD/Include)

add_library(omc::3rd::suitesparse::colamd ALIAS colamd)
target_include_directories(colamd INTERFACE SuiteSparse-5.8.1/COLAMD/Include)

add_library(omc::3rd::suitesparse::btf ALIAS btf)
target_include_directories(btf INTERFACE SuiteSparse-5.8.1/BTF/Include)

add_library(omc::3rd::suitesparse::config ALIAS suitesparseconfig)
target_include_directories(suitesparseconfig INTERFACE SuiteSparse-5.8.1/SuiteSparse_config/)


# sundials
set(KLU_LIBRARY klu)
set(AMD_LIBRARY amd)
set(COLAMD_LIBRARY colamd)
set(BTF_LIBRARY btf)
set(SUITESPARSECONFIG_LIBRARY suitesparseconfig)

option(SUNDIALS_BUILD_STATIC_LIBS "Build static libraries" ON)
option(SUNDIALS_BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(SUNDIALS_KLU_ENABLE "Enable KLU support" ON)
option(SUNDIALS_EXAMPLES_ENABLE_C "Build SUNDIALS C examples" OFF)

if(OM_OMC_ENABLE_FORTRAN)
  option(SUNDIALS_LAPACK_ENABLE "Enable Lapack support" ON)
endif()

omc_add_subdirectory(sundials-5.4.0)

## Sundials thoughtfully has organized its headers cleanly in one include/ directory
## Take advantage of that to transitively provide the headers when an external target links to
## any one of the sundials' libs.
add_library(sundials_interface INTERFACE)
target_include_directories(sundials_interface INTERFACE ${sundials_SOURCE_DIR}/include/)
## The sundials_config.h files are generated in the build directory. Add it as an include dir.
target_include_directories(sundials_interface INTERFACE ${sundials_BINARY_DIR}/include/)

## Add an interface lib for linking ot the static libs. This will transitively add
## DLINK_SUNDIALS_STATIC to anything that links to the static sundials libs.
add_library(sundials_interface_static INTERFACE)
target_link_libraries(sundials_interface_static INTERFACE sundials_interface)
target_compile_definitions(sundials_interface_static INTERFACE LINK_SUNDIALS_STATIC)

## Now that the includes and defines are attached to a utility interface library (sundials_interface_static) link it
## to the sundials static libs so they can be found when the sundials lib is linked-to from an external lib.
### Note! It should have been enough for the scope here to be INTERFACE instead of PUBLIC. However,
### that does not seem to work. This should be fine anyway.
target_link_libraries(sundials_cvode_static PUBLIC sundials_interface_static)
target_link_libraries(sundials_idas_static PUBLIC sundials_interface_static)
target_link_libraries(sundials_kinsol_static PUBLIC sundials_interface_static)
target_link_libraries(sundials_sunlinsolklu_static PUBLIC sundials_interface_static)
if(OM_OMC_ENABLE_FORTRAN)
  target_link_libraries(sundials_sunlinsollapackdense_static PUBLIC sundials_interface_static)
endif()

## Add aliases for the static libs. For readability.
add_library(omc::3rd::sundials::cvode ALIAS sundials_cvode_static)
add_library(omc::3rd::sundials::idas ALIAS sundials_idas_static)
add_library(omc::3rd::sundials::kinsol ALIAS sundials_kinsol_static)
add_library(omc::3rd::sundials::sunlinsolklu ALIAS sundials_sunlinsolklu_static)
if(OM_OMC_ENABLE_FORTRAN)
  add_library(omc::3rd::sundials::sunlinsollapackdense ALIAS sundials_sunlinsollapackdense_static)
endif()




# Ipopt
if(WITH_IPOPT)
omc_add_subdirectory(Ipopt-3.13.4)
## Just like FMIL, Ipopt assumes it will be always installed before use. So it does not organize
## its header files properly. We deal with that here.
## We collect all public headers to a new include dir in the build directory of Ipopt.

set(IPOPT_LINALG_HDRS ${IpOpt_SOURCE_DIR}/src/LinAlg/IpMatrix.hpp ${IpOpt_SOURCE_DIR}/src/LinAlg/IpSymMatrix.hpp
        ${IpOpt_SOURCE_DIR}/src/LinAlg/IpExpansionMatrix.hpp ${IpOpt_SOURCE_DIR}/src/LinAlg/IpVector.hpp
        ${IpOpt_SOURCE_DIR}/src/LinAlg/IpDenseVector.hpp ${IpOpt_SOURCE_DIR}/src/LinAlg/IpCompoundVector.hpp
        ${IpOpt_SOURCE_DIR}/src/LinAlg/IpCompoundMatrix.hpp ${IpOpt_SOURCE_DIR}/src/LinAlg/IpCompoundSymMatrix.hpp
        ${IpOpt_SOURCE_DIR}/src/LinAlg/IpSumSymMatrix.hpp ${IpOpt_SOURCE_DIR}/src/LinAlg/IpDiagMatrix.hpp
        ${IpOpt_SOURCE_DIR}/src/LinAlg/IpIdentityMatrix.hpp ${IpOpt_SOURCE_DIR}/src/LinAlg/IpScaledMatrix.hpp
        ${IpOpt_SOURCE_DIR}/src/LinAlg/IpSymScaledMatrix.hpp ${IpOpt_SOURCE_DIR}/src/LinAlg/IpZeroSymMatrix.hpp
        ${IpOpt_SOURCE_DIR}/src/LinAlg/IpBlas.hpp ${IpOpt_SOURCE_DIR}/src/LinAlg/IpLapack.hpp)

set(IPOPT_TMATRICES_HDRS ${IpOpt_SOURCE_DIR}/src/LinAlg/TMatrices/IpGenTMatrix.hpp
        ${IpOpt_SOURCE_DIR}/src/LinAlg/TMatrices/IpSymTMatrix.hpp ${IpOpt_SOURCE_DIR}/src/LinAlg/TMatrices/IpTripletHelper.hpp)

set(IPOPT_INTERFACES_HDRS ${IpOpt_SOURCE_DIR}/src/Interfaces/IpAlgTypes.hpp
        ${IpOpt_SOURCE_DIR}/src/Interfaces/IpIpoptApplication.hpp ${IpOpt_SOURCE_DIR}/src/Interfaces/IpNLP.hpp
        ${IpOpt_SOURCE_DIR}/src/Interfaces/IpReturnCodes.h ${IpOpt_SOURCE_DIR}/src/Interfaces/IpReturnCodes.hpp
        ${IpOpt_SOURCE_DIR}/src/Interfaces/IpReturnCodes_inc.h ${IpOpt_SOURCE_DIR}/src/Interfaces/IpReturnCodes.inc
        ${IpOpt_SOURCE_DIR}/src/Interfaces/IpSolveStatistics.hpp ${IpOpt_SOURCE_DIR}/src/Interfaces/IpStdCInterface.h
        ${IpOpt_SOURCE_DIR}/src/Interfaces/IpTNLP.hpp ${IpOpt_SOURCE_DIR}/src/Interfaces/IpTNLPAdapter.hpp
        ${IpOpt_SOURCE_DIR}/src/Interfaces/IpTNLPReducer.hpp)

set(IPOPT_COMMON_HDRS ${IpOpt_SOURCE_DIR}/src/Common/IpCachedResults.hpp
        ${IpOpt_SOURCE_DIR}/src/Common/IpDebug.hpp ${IpOpt_SOURCE_DIR}/src/Common/IpException.hpp
        ${IpOpt_SOURCE_DIR}/src/Common/IpJournalist.hpp ${IpOpt_SOURCE_DIR}/src/Common/IpObserver.hpp
        ${IpOpt_SOURCE_DIR}/src/Common/IpOptionsList.hpp ${IpOpt_SOURCE_DIR}/src/Common/IpoptConfig.h
        ${IpOpt_SOURCE_DIR}/src/Common/config_ipopt_default.h ${IpOpt_SOURCE_DIR}/src/Common/IpReferenced.hpp
        ${IpOpt_SOURCE_DIR}/src/Common/IpRegOptions.hpp ${IpOpt_SOURCE_DIR}/src/Common/IpSmartPtr.hpp
        ${IpOpt_SOURCE_DIR}/src/Common/IpTaggedObject.hpp ${IpOpt_SOURCE_DIR}/src/Common/IpTimedTask.hpp
        ${IpOpt_SOURCE_DIR}/src/Common/IpTypes.hpp ${IpOpt_SOURCE_DIR}/src/Common/IpUtils.hpp)

set(IPOPT_ALGORITHMS_HDRS ${IpOpt_SOURCE_DIR}/src/Algorithm/IpIpoptCalculatedQuantities.hpp
        ${IpOpt_SOURCE_DIR}/src/Algorithm/IpIpoptData.hpp ${IpOpt_SOURCE_DIR}/src/Algorithm/IpIteratesVector.hpp
        ${IpOpt_SOURCE_DIR}/src/Algorithm/IpTimingStatistics.hpp ${IpOpt_SOURCE_DIR}/src/Algorithm/IpIpoptNLP.hpp
        ${IpOpt_SOURCE_DIR}/src/Algorithm/IpOrigIpoptNLP.hpp ${IpOpt_SOURCE_DIR}/src/Algorithm/IpNLPScaling.hpp
        ${IpOpt_SOURCE_DIR}/src/Algorithm/IpAlgBuilder.hpp ${IpOpt_SOURCE_DIR}/src/Algorithm/IpIpoptAlg.hpp
        ${IpOpt_SOURCE_DIR}/src/Algorithm/IpAlgStrategy.hpp ${IpOpt_SOURCE_DIR}/src/Algorithm/IpSearchDirCalculator.hpp
        ${IpOpt_SOURCE_DIR}/src/Algorithm/IpLineSearch.hpp ${IpOpt_SOURCE_DIR}/src/Algorithm/IpMuUpdate.hpp
        ${IpOpt_SOURCE_DIR}/src/Algorithm/IpConvCheck.hpp ${IpOpt_SOURCE_DIR}/src/Algorithm/IpIterateInitializer.hpp
        ${IpOpt_SOURCE_DIR}/src/Algorithm/IpIterationOutput.hpp ${IpOpt_SOURCE_DIR}/src/Algorithm/IpHessianUpdater.hpp
        ${IpOpt_SOURCE_DIR}/src/Algorithm/IpEqMultCalculator.hpp ${IpOpt_SOURCE_DIR}/src/Algorithm/IpAugSystemSolver.hpp
        ${IpOpt_SOURCE_DIR}/src/Algorithm/IpPDSystemSolver.hpp)
set(IPOPT_LINEARSOLVERS_HDRS ${IpOpt_SOURCE_DIR}/src/Algorithm/LinearSolvers/IpSymLinearSolver.hpp)

set(IPOPT_ALL_HDRS ${IPOPT_LINALG_HDRS} ${IPOPT_TMATRICES_HDRS} ${IPOPT_INTERFACES_HDRS}
                   ${IPOPT_COMMON_HDRS} ${IPOPT_ALGORITHMS_HDRS} ${IPOPT_LINEARSOLVERS_HDRS})

file(MAKE_DIRECTORY ${IpOpt_BINARY_DIR}/include/Ipopt)
file(COPY ${IPOPT_ALL_HDRS} DESTINATION ${IpOpt_BINARY_DIR}/include/Ipopt)

# We give this new directory as include dir ONLY for targets that depend on Ipopt (not Ipopt itslef).
target_include_directories(ipopt INTERFACE ${IpOpt_BINARY_DIR}/include/)
add_library(omc::3rd::ipopt ALIAS ipopt)

endif() # WITH_IPOPT
