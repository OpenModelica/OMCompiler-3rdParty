# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
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
# examples/ida level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# C examples
if(SUNDIALS_EXAMPLES_ENABLE_C)
  add_subdirectory(serial)
  if(SUNDIALS_OPENMP_ENABLE AND OPENMP_FOUND)
    add_subdirectory(C_openmp)
  endif()
  if(SUNDIALS_MPI_ENABLE AND MPI_C_FOUND)
    add_subdirectory(parallel)
  endif()
  if(PETSC_ENABLE AND PETSC_FOUND)
    add_subdirectory(petsc)
  endif()
endif()

# C++ examples
if(SUNDIALS_EXAMPLES_ENABLE_CXX)
  if(Trilinos_ENABLE AND SUNDIALS_TRILINOS)
    add_subdirectory(trilinos)
  endif()
endif()

# Fortran examples
if(F77_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F77 AND F77_FOUND)
  add_subdirectory(fcmix_serial)
  if(SUNDIALS_OPENMP_ENABLE AND OPENMP_FOUND)
    add_subdirectory(fcmix_openmp)
  endif()
  if(SUNDIALS_PTHREAD_ENABLE AND PTHREADS_FOUND)
    add_subdirectory(fcmix_pthreads)
  endif()
  if(SUNDIALS_MPI_ENABLE AND MPI_Fortran_FOUND)
    add_subdirectory(fcmix_parallel)
  endif()
endif()

if(F2003_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F2003 AND F2003_FOUND)
  add_subdirectory(F2003_serial)
endif()

# cuda examples
if(EXAMPLES_ENABLE_CUDA)

  if(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
    add_subdirectory(cuda)
    if(SUNDIALS_MPI_ENABLE AND MPI_C_FOUND)
      add_subdirectory(mpicuda)
    endif()
  endif()

  if(RAJA_ENABLE AND RAJA_FOUND)
    add_subdirectory(raja)
    if(SUNDIALS_MPI_ENABLE AND MPI_C_FOUND)
      add_subdirectory(mpiraja)
    endif()
  endif()

endif(EXAMPLES_ENABLE_CUDA)
