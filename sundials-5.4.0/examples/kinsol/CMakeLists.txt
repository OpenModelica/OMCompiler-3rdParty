# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
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
# examples/kinsol level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# C examples
if(SUNDIALS_EXAMPLES_ENABLE_C)
  add_subdirectory(serial)
  if(SUNDIALS_OPENMP_ENABLE AND OPENMP_FOUND)
    # the only example here need special handling from testrunner (not yet implemented)
    add_subdirectory(C_openmp)
  endif()
  if(SUNDIALS_MPI_ENABLE AND MPI_C_FOUND)
    add_subdirectory(parallel)
  endif()
endif()

# Fortran examples
if(F77_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F77 AND F77_FOUND)
  add_subdirectory(fcmix_serial)
  if(SUNDIALS_MPI_ENABLE AND MPI_Fortran_FOUND)
    add_subdirectory(fcmix_parallel)
  endif()
endif()

if(F2003_INTERFACE_ENABLE AND EXAMPLES_ENABLE_F2003 AND F2003_FOUND)
  add_subdirectory(F2003_serial)
endif()
