project(OMCompiler_3rdParty)

# antlr
omc_add_subdirectory(antlr)
message(STATUS ${OMAntlr3_ANTLRJAR})

# lpsolve55
omc_add_subdirectory(lpsolve)


# Boehm GC
omc_add_subdirectory(gc)

# Intel TBB
set(TBB_BUILD_TESTS OFF CACHE BOOL "Build TBB tests and enable testing infrastructure")
set(TBB_BUILD_SHARED OFF CACHE BOOL "Build TBB shared library")
set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL "Build TBB malloc library")
set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL "Build TBB malloc proxy library")
omc_add_subdirectory(tbb_cmake)

# sundials
omc_add_subdirectory(sundials)

# CMinpack
omc_add_subdirectory(CMinpack)

# # cppzmq
# omc_add_subdirectory(cppzmq)

# FMIL
set (FMILIB_GENERATE_DOXYGEN_DOC OFF CACHE BOOL "Generate doxygen doc target")
set (FMILIB_BUILD_TESTS OFF CACHE BOOL "Build tests")
omc_add_subdirectory(FMIL)
omc_add_to_report(FMILIB_FMI_PLATFORM)


# metis
omc_add_subdirectory(metis-5.1.0)

# # regex
# omc_add_subdirectory(regex-0.12)

# SuiteSparse
omc_add_subdirectory(SuiteSparse)
