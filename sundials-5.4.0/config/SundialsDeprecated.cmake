# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Print warning is the user sets a deprecated CMake variable and
# copy the value into the correct CMake variable
# ---------------------------------------------------------------

# macro to print warning for deprecated CMake variable
macro(PRINT_DEPRECATED old_variable new_variable)
  print_warning("${old_variable} is deprecated and will be removed in the future."
                "Copying value to ${new_variable}.")
endmacro()

if(DEFINED SUNDIALS_EXAMPLES_ENABLE)
  print_deprecated(SUNDIALS_EXAMPLES_ENABLE SUNDIALS_EXAMPLES_ENABLE_C)
  force_variable(SUNDIALS_EXAMPLES_ENABLE_C BOOL "Build SUNDIALS C examples" ${SUNDIALS_EXAMPLES_ENABLE})
  unset(SUNDIALS_EXAMPLES_ENABLE CACHE)
endif()

if(DEFINED CXX_ENABLE)
  print_deprecated(CXX_ENABLE SUNDIALS_EXAMPLES_ENABLE_CXX)
  force_variable(SUNDIALS_EXAMPLES_ENABLE_CXX BOOL "Build ARKode C++ examples" ${CXX_ENABLE})
  unset(CXX_ENABLE CACHE)
endif()

if(DEFINED F90_ENABLE)
  print_deprecated(F90_ENABLE EXAMPLES_ENABLE_F90)
  force_variable(EXAMPLES_ENABLE_F90 BOOL "Build ARKode Fortran90 examples" ${F90_ENABLE})
  unset(F90_ENABLE CACHE)
endif()

if(DEFINED FCMIX_ENABLE)
  print_deprecated(FCMIX_ENABLE F77_INTERFACE_ENABLE)
  force_variable(F77_INTERFACE_ENABLE BOOL "Build Fortran 77 interfaces" ${FCMIX_ENABLE})
  unset(FCMIX_ENABLE CACHE)
endif()

# SUNDIALS_INDEX_TYPE got new behavior
if(SUNDIALS_INDEX_TYPE)
  string(TOUPPER ${SUNDIALS_INDEX_TYPE} tmp)

  if(tmp STREQUAL "INT32_T")
    print_warning("SUNDIALS_INDEX_TYPE overrides the standard types SUNDIALS looks for."
    "Setting SUNDIALS_INDEX_SIZE to 32 and clearing SUNDIALS_INDEX_TYPE.")
    force_variable(SUNDIALS_INDEX_SIZE STRING "SUNDIALS index size" 32)
    force_variable(SUNDIALS_INDEX_TYPE STRING "SUNDIALS index type" "")
  elseif(tmp STREQUAL "INT64_T")
    print_warning("SUNDIALS_INDEX_TYPE overrides the standard types SUNDIALS looks for."
    "Setting SUNDIALS_INDEX_SIZE to 64 and clearing SUNDIALS_INDEX_TYPE.")
    force_variable(SUNDIALS_INDEX_SIZE STRING "SUNDIALS index size" 64)
    force_variable(SUNDIALS_INDEX_TYPE STRING "SUNDIALS index type" "")
  else()
    print_warning("SUNDIALS_INDEX_TYPE overrides the standard types SUNDIALS looks for." "")
  endif()
endif()

if(DEFINED MPI_MPICC)
  print_deprecated(MPI_MPICC MPI_C_COMPILER)
  force_variable(MPI_C_COMPILER FILEPATH "MPI C compiler" ${MPI_MPICC})
  unset(MPI_MPICC CACHE)
endif()

if(DEFINED MPI_MPICXX)
  print_deprecated(MPI_MPICXX MPI_CXX_COMPILER)
  force_variable(MPI_CXX_COMPILER FILPATH "MPI C++ compiler" ${MPI_MPICXX})
  unset(MPI_MPICXX CACHE)
endif()

if((DEFINED MPI_MPIF77) OR (DEFINED MPI_MPIF90))
  if(DEFINED MPI_MPIF90)
    print_warning("MPI_MPIF77 and MPI_MPIF90 are deprecated and will be removed in the future." "Copying MPI_MPIF90 value to MPI_Fortran_COMPILER")
    force_variable(MPI_Fortran_COMPILER FILEPATH "MPI Fortran compiler" ${MPI_MPIF90})
  else()
    print_warning("MPI_MPIF77 and MPI_MPIF90 are deprecated and will be removed in the future." "Copying MPI_MPIF77 value to MPI_Fortran_COMPILER")
    force_variable(MPI_Fortran_COMPILER FILEPATH "MPI Fortran compiler" ${MPI_MPIF77})
  endif()
  unset(MPI_MPIF77 CACHE)
  unset(MPI_MPIF90 CACHE)
endif()

if(DEFINED MPI_RUN_COMMAND)
  print_deprecated(MPI_RUN_COMMAND MPIEXEC_EXECUTABLE)
  force_variable(MPIEXEC_EXECUTABLE FILEPATH "MPI run command" ${MPI_RUN_COMMAND})
  unset(MPI_RUN_COMMAND CACHE)
endif()


###############################################################################
# Secret option to install impl header files.                                 #
# TODO: remove after Sept. 2019                                               #
###############################################################################
if(DEFINED _INSTALL_IMPL_FILES)
  if(BUILD_ARKODE)
    install(FILES ${PROJECT_SOURCE_DIR}/src/arkode/arkode_impl.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arkode)
  endif()
  if(BUILD_CVODE)
    install(FILES ${PROJECT_SOURCE_DIR}/src/cvode/cvode_impl.h   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cvode)
  endif()
  if(BUILD_CVODES)
    install(FILES ${PROJECT_SOURCE_DIR}/src/cvodes/cvodes_impl.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cvodes)
  endif()
  if(BUILD_IDA)
    install(FILES ${PROJECT_SOURCE_DIR}/src/ida/ida_impl.h       DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ida)
  endif()
  if(BUILD_IDAS)
    install(FILES ${PROJECT_SOURCE_DIR}/src/idas/idas_impl.h     DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/idas)
  endif()
  if(BUILD_KINSOL)
    install(FILES ${PROJECT_SOURCE_DIR}/src/kinsol/kinsol_impl.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/kinsol)
  endif()
endif()
