# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# linear solver level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# Always add SUNDIALS provided linear solver modules
add_subdirectory(band)
add_subdirectory(dense)
add_subdirectory(pcg)

add_subdirectory(spbcgs)
add_subdirectory(spfgmr)
add_subdirectory(spgmr)
add_subdirectory(sptfqmr)

# Add linear solver modules with TPL dependencies
if(SUNDIALS_KLU_ENABLE AND KLU_FOUND)
  add_subdirectory(klu)
endif(SUNDIALS_KLU_ENABLE AND KLU_FOUND)

if(SUNDIALS_SUPERLUDIST)
  add_subdirectory(superludist)
endif()

if(SUNDIALS_PTHREAD_ENABLE AND SUPERLUMT_FOUND)
  add_subdirectory(superlumt)
endif(SUNDIALS_PTHREAD_ENABLE AND SUPERLUMT_FOUND)

if(SUNDIALS_LAPACK_ENABLE AND LAPACK_FOUND)
  add_subdirectory(lapackband)
  add_subdirectory(lapackdense)
endif(SUNDIALS_LAPACK_ENABLE AND LAPACK_FOUND)

if(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
  if(SUNDIALS_INDEX_SIZE MATCHES "32")
    add_subdirectory(cusolversp)
  endif()
endif(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
