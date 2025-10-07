# This file is copied into MUMPS folder once the library downloaded
#

cmake_minimum_required (VERSION 3.10)

project (MUMPS C CXX Fortran)

set(MUMPS_LAPACK_LIB_PATH "None" CACHE PATH "The Lapack library library Path")

message(STATUS "Building MUMPS with METIS: ${MUMPS_USE_METIS}")

set(INCLUDEDIR ${CMAKE_CURRENT_BINARY_DIR}/MUMPS/include)

#------------------------------------------------------------
# Detect 64 bits
#------------------------------------------------------------

if (CMAKE_SIZEOF_VOID_P EQUAL 4)
  set(HAVE_64_BIT 0)
else ()
  set(HAVE_64_BIT 1)
endif ()

#------------------------------------------------------------
# MPI
#------------------------------------------------------------

if (NOT MUMPS_USE_LIBSEQ AND NOT MPI_FOUND)
  # Use MPI_EXTRA_LIBNAMES to specify required libraries under MPI
  set(MPI_EXTRA_LIBNAMES fmpich mpich)
  find_package(MPI)
endif ()

if (NOT MUMPS_USE_LIBSEQ)
  add_definitions("-Dpord")
endif ()

#------------------------------------------------------------
# ATLAS / MKL
#------------------------------------------------------------

# Use MKL_LIBS to specify the required libs under MKL
if (WIN32)
  if (HAVE_64_BIT)
    set(MKL_LIBS mkl_lapack95_lp64
                 mkl_blas95_lp64
                 mkl_solver_lp64
                 mkl_intel_lp64
                 mkl_intel_thread
                 mkl_core
                 mkl_blacs_lp64_dll)
  else ()
    set(MKL_LIBS mkl_lapack95
                 mkl_blas95
                 mkl_solver
                 mkl_intel_c
                 mkl_intel_thread
                 mkl_core
                 mkl_blacs_dll)
  endif ()
else ()
  set(MKL_LIBS mkl_lapack95_lp64
               mkl_blas95_lp64
               mkl_solver_lp64
               mkl_intel_lp64
               mkl_intel_thread
               mkl_core
               mkl_blacs_openmpi_lp64)
endif ()


if (NOT LAPACK_FOUND AND NOT DOWNLOAD_LAPACK)
  message(STATUS "Searching for LAPACK")
  find_package(LAPACK REQUIRED)
endif ()

#------------------------------------------------------------
# configure config*.h
#------------------------------------------------------------

# if MPI found, then do not compile libseq ...

if (IPOPT_INT64)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/mumps_int_def64_h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/mumps_int_def.h)
else ()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/mumps_int_def32_h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/mumps_int_def.h)
endif ()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif ()

if (MUMPS_USE_LIBSEQ)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libseq)
endif ()

if (MUMPS_USE_METIS)
  include_directories(${METIS_INCLUDE_DIR})
endif ()

if (MKL_FOUND)
  link_directories(${MKL_PATH})
endif ()

if (MPI_FOUND)
  link_directories(${MPI_LIBRARY_PATH})
endif ()

#------------------------------------------------------------
# METIS
#------------------------------------------------------------

if (MUMPS_USE_METIS)
  add_definitions(-Dmetis)
endif ()

if (NOT MUMPS_LAPACK_LIB_PATH STREQUAL "None")
  link_directories(${MUMPS_LAPACK_LIB_PATH})
endif ()

#------------------------------------------------------------
# Build
#------------------------------------------------------------

set(MUMPS_PORD_SRCS PORD/lib/graph.c
                    PORD/lib/gbipart.c
                    PORD/lib/gbisect.c
                    PORD/lib/ddcreate.c
                    PORD/lib/ddbisect.c
                    PORD/lib/nestdiss.c
                    PORD/lib/multisector.c
                    PORD/lib/gelim.c
                    PORD/lib/bucket.c
                    PORD/lib/tree.c
                    PORD/lib/symbfac.c
                    PORD/lib/interface.c
                    PORD/lib/sort.c
                    PORD/lib/minpriority.c)

set(MUMPS_LIBSEQ_SRCS libseq/mpi.f
                      libseq/mpic.c
                      libseq/elapse.c)

set(MUMPS_COMMON_SRCS src/ana_AMDMF.F
                  src/ana_blk.F
                  src/ana_blk_m.F
                  src/ana_omp_m.F
                  src/ana_orderings.F
                  src/ana_orderings_wrappers_m.F
                  src/ana_set_ordering.F
                  src/bcast_errors.F
                  src/double_linked_list.F
                  src/estim_flops.F
                  src/fac_asm_build_sort_index_ELT_m.F
                  src/fac_asm_build_sort_index_m.F
                  src/fac_descband_data_m.F
                  src/fac_future_niv2_mod.F
                  src/fac_maprow_data_m.F
                  src/front_data_mgt_m.F
                  src/lr_common.F
                  src/lr_stats.F
                  src/mumps_addr.c
                  src/mumps_comm_buffer_common.F
                  src/mumps_common.c
                  src/mumps_config_file_C.c
                  src/mumps_flytes.c
                  src/mumps_intr_types_common.F
                  src/mumps_io_basic.c
                  src/mumps_io.c
                  src/mumps_io_err.c
                  src/mumps_io_thread.c
                  src/mumps_l0_omp_m.F
                  src/mumps_load.F
                  src/mumps_memory_mod.F
                  src/mumps_metis64.c
                  src/mumps_metis.c
                  src/mumps_metis_int.c
                  src/mumps_mpitoomp_m.F
                  src/mumps_numa.c
                  src/mumps_ooc_common.F
                  src/mumps_pivnul_mod.F
                  src/mumps_pord.c
                  src/mumps_print_defined.F
                  src/mumps_save_restore_C.c
                  src/mumps_scotch64.c
                  src/mumps_scotch.c
                  src/mumps_scotch_int.c
                  src/mumps_static_mapping.F
                  src/mumps_thread_affinity.c
                  src/mumps_register_thread.c
                  src/mumps_thread.c
                  src/mumps_type2_blocking.F
                  src/mumps_version.F
                  src/omp_tps_common_m.F
                  src/sol_common.F
                  src/sol_ds_common_m.F
                  src/sol_omp_common_m.F
                  src/tools_common_m.F
                  src/tools_common.F
                  libseq/mpi.f
                  libseq/mpic.c)

set(MUMPS_C_SRCS  src/mumps_c.c
                  src/cana_aux_ELT.F
                  src/cana_aux.F
                  src/cana_aux_par.F
                  src/cana_dist_m.F
                  src/cana_driver.F
                  src/cana_LDLT_preprocess.F
                  src/cana_lr.F
                  src/cana_mtrans.F
                  src/cana_reordertree.F
                  src/carrowheads.F
                  src/cbcast_int.F
                  src/cend_driver.F
                  src/cfac_asm_ELT.F
                  src/cfac_asm.F
                  src/cfac_asm_master_ELT_m.F
                  src/cfac_asm_master_m.F
                  src/cfac_compact_factors_m.F
                  src/cfac_b.F
                  src/cfac_diag.F
                  src/cfac_determinant.F
                  src/cfac_distrib_distentry.F
                  src/cfac_distrib_ELT.F
                  src/cfac_dist_arrowheads_omp.F
                  src/cfac_driver.F
                  src/cfac_front_aux.F
                  src/cfac_front_LU_type1.F
                  src/cfac_front_LU_type2.F
                  src/cfac_front_type2_aux.F
                  src/cfac_lastrtnelind.F
                  src/cfac_lr.F
                  src/cfac_mem_alloc_cb.F
                  src/cfac_mem_compress_cb.F
                  src/cfac_mem_dynamic.F
                  src/cfac_mem_free_block_cb.F
                  src/cfac_mem_stack_aux.F
                  src/cfac_mem_stack.F
                  src/cfac_omp_m.F
                  src/cfac_front_LDLT_type1.F
                  src/cfac_front_LDLT_type2.F
                  src/cfac_par_m.F
                  src/cfac_process_band.F
                  src/cfac_process_bf.F
                  src/cfac_process_blfac_slave.F
                  src/cfac_process_blocfacto.F
                  src/cfac_process_contrib_type1.F
                  src/cfac_process_contrib_type2.F
                  src/cfac_process_contrib_type3.F
                  src/cfac_process_end_facto_slave.F
                  src/cfac_process_maprow.F
                  src/cfac_process_master2.F
                  src/cfac_process_message.F
                  src/cfac_process_root2slave.F
                  src/cfac_process_root2son.F
                  src/cfac_process_rtnelind.F
                  src/cfac_root_parallel.F
                  src/cfac_scalings.F
                  src/cfac_determinant.F
                  src/cfac_scalings_simScaleAbs.F
                  src/cfac_scalings_simScale_util.F
                  src/cfac_sol_pool.F
                  src/cfac_type3_symmetrize.F
                  src/cfac_sol_l0omp_m.F
                  src/cfac_sispointers_m.F
                  src/cini_defaults.F
                  src/cini_driver.F
                  src/clr_core.F
                  src/clr_type.F
                  src/cmumps_comm_buffer.F
                  src/cmumps_config_file.F
                  src/cmumps_driver.F
                  src/cmumps_f77.F
                  src/cmumps_gpu.c
                  src/cmumps_intr_types.F
                  src/cmumps_iXamax.F
                  src/cmumps_lr_data_m.F
                  src/cmumps_mpi3_mod.F
                  src/cmumps_ooc_buffer.F
                  src/cmumps_ooc.F
                  src/cmumps_save_restore.F
                  src/cmumps_save_restore_files.F
                  src/cmumps_sol_es.F
                  src/cmumps_struc_def.F
                  src/comp_tps_m.F
                  src/cooc_panel_piv.F
                  src/crank_revealing.F
                  src/csol_aux.F
                  src/csol_bwd_aux.F
                  src/csol_bwd.F
                  src/csol_c.F
                  src/csol_distrhs.F
                  src/csol_distsol.F
                  src/csol_driver.F
                  src/csol_fwd_aux.F
                  src/csol_fwd.F
                  src/csol_lr.F
                  src/csol_matvec.F
                  src/csol_omp_m.F
                  src/csol_root_parallel.F
                  src/cstatic_ptr_m.F
                  src/ctools.F
                  src/ctype3_root.F)

set(MUMPS_D_SRCS  src/mumps_c.c
                  src/dana_aux_ELT.F
                  src/dana_aux.F
                  src/dana_aux_par.F
                  src/dana_dist_m.F
                  src/dana_driver.F
                  src/dana_LDLT_preprocess.F
                  src/dana_lr.F
                  src/dana_mtrans.F
                  src/dana_reordertree.F
                  src/darrowheads.F
                  src/dbcast_int.F
                  src/dend_driver.F
                  src/dfac_asm_ELT.F
                  src/dfac_asm.F
                  src/dfac_asm_master_ELT_m.F
                  src/dfac_asm_master_m.F
                  src/dfac_compact_factors_m.F
                  src/dfac_b.F
                  src/dfac_diag.F
                  src/dfac_determinant.F
                  src/dfac_distrib_distentry.F
                  src/dfac_distrib_ELT.F
                  src/dfac_dist_arrowheads_omp.F
                  src/dfac_driver.F
                  src/dfac_front_aux.F
                  src/dfac_front_LDLT_type1.F
                  src/dfac_front_LDLT_type2.F
                  src/dfac_front_LU_type1.F
                  src/dfac_front_LU_type2.F
                  src/dfac_front_type2_aux.F
                  src/dfac_lastrtnelind.F
                  src/dfac_lr.F
                  src/dfac_mem_alloc_cb.F
                  src/dfac_mem_compress_cb.F
                  src/dfac_mem_dynamic.F
                  src/dfac_mem_free_block_cb.F
                  src/dfac_mem_stack_aux.F
                  src/dfac_mem_stack.F
                  src/dfac_omp_m.F
                  src/dfac_par_m.F
                  src/dfac_process_band.F
                  src/dfac_process_bf.F
                  src/dfac_process_blfac_slave.F
                  src/dfac_process_blocfacto.F
                  src/dfac_process_blocfacto_LDLT.F
                  src/dfac_process_contrib_type1.F
                  src/dfac_process_contrib_type2.F
                  src/dfac_process_contrib_type3.F
                  src/dfac_process_end_facto_slave.F
                  src/dfac_process_maprow.F
                  src/dfac_process_master2.F
                  src/dfac_process_message.F
                  src/dfac_process_root2slave.F
                  src/dfac_process_root2son.F
                  src/dfac_process_rtnelind.F
                  src/dfac_root_parallel.F
                  src/dfac_scalings.F
                  src/dfac_determinant.F
                  src/dfac_scalings_simScaleAbs.F
                  src/dfac_scalings_simScale_util.F
                  src/dfac_sol_l0omp_m.F
                  src/dfac_sispointers_m.F
                  src/dfac_sol_pool.F
                  src/dfac_type3_symmetrize.F
                  src/dini_defaults.F
                  src/dini_driver.F
                  src/dlr_core.F
                  src/dlr_type.F
                  src/dmumps_comm_buffer.F
                  src/dmumps_config_file.F
                  src/dmumps_driver.F
                  src/dmumps_f77.F
                  src/dmumps_gpu.c
                  src/dmumps_intr_types.F
                  src/dmumps_iXamax.F
                  src/dmumps_lr_data_m.F
                  src/dmumps_mpi3_mod.F
                  src/dmumps_ooc_buffer.F
                  src/dmumps_ooc.F
                  src/dmumps_save_restore.F
                  src/dmumps_save_restore_files.F
                  src/dmumps_sol_es.F
                  src/dmumps_struc_def.F
                  src/domp_tps_m.F
                  src/dooc_panel_piv.F
                  src/drank_revealing.F
                  src/dsol_aux.F
                  src/dsol_bwd_aux.F
                  src/dsol_bwd.F
                  src/dsol_c.F
                  src/dsol_distrhs.F
                  src/dsol_distsol.F
                  src/dsol_driver.F
                  src/dsol_fwd_aux.F
                  src/dsol_fwd.F
                  src/dsol_lr.F
                  src/dsol_matvec.F
                  src/dsol_omp_m.F
                  src/dsol_root_parallel.F
                  src/dstatic_ptr_m.F
                  src/dtools.F
                  src/dtype3_root.F)

set(MUMPS_S_SRCS src/mumps_c.c
                src/sana_aux_ELT.F
                src/sana_aux.F
                src/sana_aux_par.F
                src/sana_dist_m.F
                src/sana_driver.F
                src/sana_LDLT_preprocess.F
                src/sana_lr.F
                src/sana_mtrans.F
                src/sana_reordertree.F
                src/sarrowheads.F
                src/sbcast_int.F
                src/send_driver.F
                src/sfac_asm_ELT.F
                src/sfac_asm.F
                src/sfac_asm_master_ELT_m.F
                src/sfac_asm_master_m.F
                src/sfac_compact_factors_m.F
                src/sfac_b.F
                src/sfac_diag.F
                src/sfac_determinant.F
                src/sfac_distrib_distentry.F
                src/sfac_distrib_ELT.F
                src/sfac_dist_arrowheads_omp.F
                src/sfac_driver.F
                src/sfac_front_aux.F
                src/sfac_front_LDLT_type1.F
                src/sfac_front_LDLT_type2.F
                src/sfac_front_LU_type1.F
                src/sfac_front_LU_type2.F
                src/sfac_front_type2_aux.F
                src/sfac_lastrtnelind.F
                src/sfac_lr.F
                src/sfac_mem_alloc_cb.F
                src/sfac_mem_compress_cb.F
                src/sfac_mem_dynamic.F
                src/sfac_mem_free_block_cb.F
                src/sfac_mem_stack_aux.F
                src/sfac_mem_stack.F
                src/sfac_omp_m.F
                src/sfac_par_m.F
                src/sfac_process_band.F
                src/sfac_process_bf.F
                src/sfac_process_blfac_slave.F
                src/sfac_process_blocfacto.F
                src/sfac_process_blocfacto_LDLT.F
                src/sfac_process_contrib_type1.F
                src/sfac_process_contrib_type2.F
                src/sfac_process_contrib_type3.F
                src/sfac_process_end_facto_slave.F
                src/sfac_process_maprow.F
                src/sfac_process_master2.F
                src/sfac_process_message.F
                src/sfac_process_root2slave.F
                src/sfac_process_root2son.F
                src/sfac_process_rtnelind.F
                src/sfac_root_parallel.F
                src/sfac_scalings.F
                src/sfac_determinant.F
                src/sfac_scalings_simScaleAbs.F
                src/sfac_scalings_simScale_util.F
                src/sfac_sol_pool.F
                src/sfac_type3_symmetrize.F
                src/sini_defaults.F
                src/sini_driver.F
                src/slr_core.F
                src/slr_type.F
                src/smumps_comm_buffer.F
                src/smumps_config_file.F
                src/smumps_driver.F
                src/smumps_f77.F
                src/smumps_gpu.c
                src/smumps_intr_types.F
                src/smumps_iXamax.F
                src/smumps_lr_data_m.F
                src/smumps_mpi3_mod.F
                src/smumps_ooc_buffer.F
                src/smumps_ooc.F
                src/smumps_save_restore.F
                src/smumps_save_restore_files.F
                src/smumps_sol_es.F
                src/smumps_struc_def.F
                src/somp_tps_m.F
                src/sooc_panel_piv.F
                src/srank_revealing.F
                src/ssol_aux.F
                src/ssol_bwd_aux.F
                src/ssol_bwd.F
                src/ssol_c.F
                src/ssol_distrhs.F
                src/ssol_distsol.F
                src/ssol_driver.F
                src/ssol_fwd_aux.F
                src/ssol_fwd.F
                src/ssol_lr.F
                src/ssol_matvec.F
                src/ssol_omp_m.F
                src/ssol_root_parallel.F
                src/sstatic_ptr_m.F
                src/stools.F
                src/stype3_root.F)

set(MUMPS_Z_SRCS src/mumps_c.c
                src/zana_aux_ELT.F
                src/zana_aux.F
                src/zana_aux_par.F
                src/zana_dist_m.F
                src/zana_driver.F
                src/zana_LDLT_preprocess.F
                src/zana_lr.F
                src/zana_mtrans.F
                src/zana_reordertree.F
                src/zarrowheads.F
                src/zbcast_int.F
                src/zend_driver.F
                src/zfac_asm_ELT.F
                src/zfac_asm.F
                src/zfac_asm_master_ELT_m.F
                src/zfac_asm_master_m.F
                src/zfac_compact_factors_m.F
                src/zfac_b.F
                src/zfac_diag.F
                src/zfac_determinant.F
                src/zfac_distrib_distentry.F
                src/zfac_distrib_ELT.F
                src/zfac_dist_arrowheads_omp.F
                src/zfac_driver.F
                src/zfac_front_aux.F
                src/zfac_front_LU_type1.F
                src/zfac_front_LU_type2.F
                src/zfac_front_type2_aux.F
                src/zfac_lastrtnelind.F
                src/zfac_lr.F
                src/zfac_mem_alloc_cb.F
                src/zfac_mem_compress_cb.F
                src/zfac_mem_dynamic.F
                src/zfac_mem_free_block_cb.F
                src/zfac_mem_stack_aux.F
                src/zfac_mem_stack.F
                src/zfac_omp_m.F
                src/zfac_par_m.F
                src/zfac_process_band.F
                src/zfac_process_bf.F
                src/zfac_process_blfac_slave.F
                src/zfac_process_blocfacto.F
                src/zfac_process_blocfacto_LDLT.F
                src/zfac_process_contrib_type1.F
                src/zfac_process_contrib_type2.F
                src/zfac_process_contrib_type3.F
                src/zfac_process_end_facto_slave.F
                src/zfac_process_maprow.F
                src/zfac_process_master2.F
                src/zfac_process_message.F
                src/zfac_process_root2slave.F
                src/zfac_process_root2son.F
                src/zfac_process_rtnelind.F
                src/zfac_root_parallel.F
                src/zfac_scalings.F
                src/zfac_determinant.F
                src/zfac_scalings_simScaleAbs.F
                src/zfac_scalings_simScale_util.F
                src/zfac_sol_pool.F
                src/zfac_type3_symmetrize.F
                src/zini_defaults.F
                src/zini_driver.F
                src/zlr_core.F
                src/zlr_type.F
                src/zmumps_comm_buffer.F
                src/zmumps_config_file.F
                src/zmumps_driver.F
                src/zmumps_f77.F
                src/zmumps_gpu.c
                src/zmumps_intr_types.F
                src/zmumps_iXamax.F
                src/zmumps_lr_data_m.F
                src/zmumps_mpi3_mod.F
                src/zmumps_ooc_buffer.F
                src/zmumps_ooc.F
                src/zmumps_save_restore.F
                src/zmumps_save_restore_files.F
                src/zmumps_sol_es.F
                src/zmumps_struc_def.F
                src/zomp_tps_m.F
                src/zooc_panel_piv.F
                src/zrank_revealing.F
                src/zsol_aux.F
                src/zsol_bwd_aux.F
                src/zsol_bwd.F
                src/zsol_c.F
                src/zsol_distrhs.F
                src/zsol_distsol.F
                src/zsol_driver.F
                src/zsol_fwd_aux.F
                src/zsol_fwd.F
                src/zsol_lr.F
                src/zsol_matvec.F
                src/zsol_omp_m.F
                src/zsol_root_parallel.F
                src/zstatic_ptr_m.F
                src/ztools.F
                src/ztype3_root.F)

set(LINK_LIBS ${LINK_LIBS} LAPACK_TARGET)

add_library(mumps_flags INTERFACE)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(mumps_flags INTERFACE
        $<$<COMPILE_LANGUAGE:Fortran>:
        -w
        -fcray-pointer
        -fallow-argument-mismatch
        -fall-intrinsics
        -finit-local-zero>
    )

    target_link_libraries(mumps_flags INTERFACE gfortran)
endif ()

target_compile_definitions(mumps_flags INTERFACE
    ALLOW_NON_INIT
    intel_
    Add_
    MUMPS_ARITH=MUMPS_ARITH_d
)

if (WIN32)
  if (MUMPS_USE_LIBSEQ)
    set(LINK_LIBS ${LINK_LIBS} seq)
  else ()
    set(LINK_LIBS ${LINK_LIBS} ${MPI_LIB_NAME})
  endif ()

  if (MKL_FOUND)
    set(LINK_LIBS ${LINK_LIBS} ${ATLAS_LIBRARIES} libiomp5mt)
  endif ()
else ()
  if (MKL_FOUND)
    set(LINK_LIBS ${LINK_LIBS} ${ATLAS_LIBRARIES} iomp5)
  endif ()

  if (MUMPS_USE_LIBSEQ)
    set(LINK_LIBS ${LINK_LIBS} seq)
  else ()
    set(LINK_LIBS ${LINK_LIBS} ${MPI_LIB_NAME})
  endif ()

  set(LINK_LIBS ${LINK_LIBS} pthread)
endif ()

if (MUMPS_USE_METIS)
  set(LINK_LIBS ${LINK_LIBS} METIS_TARGET m)
endif ()

if (MUMPS_USE_LIBSEQ)
  add_library(seq STATIC ${MUMPS_LIBSEQ_SRCS})
  target_link_libraries(seq PRIVATE std_global_flags mumps_flags)
endif ()

if (NOT MUMPS_USE_LIBSEQ)
  add_library(pord STATIC ${MUMPS_PORD_SRCS})
  target_include_directories(pord BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/PORD/include)
endif ()

# To allow the link of examples on the cluster
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(LINK_LIBS ${LINK_LIBS}
                ifcore)
endif ()

add_library(mumps_common STATIC ${MUMPS_COMMON_SRCS})
target_link_libraries(mumps_common PRIVATE std_global_flags mumps_flags)
target_link_libraries(mumps_common PUBLIC ${LINK_LIBS})

# TODO: add MUMPS single precision

# add_library(smumps STATIC ${MUMPS_S_SRCS})
# target_link_libraries(smumps mumps_common)
#

add_library(dmumps STATIC ${MUMPS_D_SRCS})
target_link_libraries(dmumps PRIVATE std_global_flags mumps_flags)
target_link_libraries(dmumps PUBLIC mumps_common ${LINK_LIBS})
target_include_directories(dmumps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(dmumps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libseq)

# add_library(cmumps STATIC ${MUMPS_C_SRCS})
# target_link_libraries(cmumps mumps_common)

# add_library(zmumps STATIC ${MUMPS_Z_SRCS})
# target_link_libraries(zmumps mumps_common)

add_executable(dsimple_test examples/dsimpletest.F)
target_link_libraries(dsimple_test dmumps mumps_common ${LINK_LIBS})
set_target_properties(dsimple_test PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(c_example examples/c_example.c)
target_link_libraries(c_example dmumps mumps_common ${LINK_LIBS})
if (WIN32)
  # Under windows, this line is required to allow compilation of the MUMPS C example.
  # Under Linux, this line makes the link phase hangs because of multiply defined main symbol
  set_target_properties(c_example PROPERTIES LINKER_LANGUAGE Fortran)
endif ()


# Install rules
if (MUMPS_USE_LIBSEQ)
  install(TARGETS seq)
endif ()

if (NOT MUMPS_USE_LIBSEQ)
  install(TARGETS pord)
endif ()

install(TARGETS mumps_common)

install(TARGETS dmumps PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mumps)

if (MUMPS_USE_LIBSEQ)
  install(DIRECTORY libseq/
          DESTINATION ${INCLUDEDIR}/
          PATTERN "*.h")

  install(FILES libseq/mpi.h
          DESTINATION ${INCLUDEDIR}/
          RENAME mumps_mpi.h)
else ()
  install(FILES ${MPI_C_INCLUDE_PATH}/libseq/mpi.h
          DESTINATION ${INCLUDEDIR}/
          RENAME mumps_mpi.h)
endif ()

install(DIRECTORY include/
        DESTINATION ${INCLUDEDIR}/
        PATTERN "*.h")

install(DIRECTORY PORD/include/
        DESTINATION ${INCLUDEDIR}/
        PATTERN "*.h")
