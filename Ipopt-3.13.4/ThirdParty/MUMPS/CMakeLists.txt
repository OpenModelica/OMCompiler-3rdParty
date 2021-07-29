
cmake_minimum_required (VERSION 2.8)
project (MUMPS C CXX Fortran)

## The OM Ipopt uses libseq always. I think libseq contains stubs for MPI.
## OM does not use MPI as far as I know so we are good. The alternative to
## libseq is maybe to use MPI. This is what I have understood so far.
# option(MUMPS_USE_LIBSEQ   "Use the MUMPS sequential MPI stub" ON)

option(MUMPS_BUILD_SHARED_LIBS    "Build libraries as shared libraries"  OFF)
option(MUMPS_USE_METIS    "Use the Metis library"             OFF)

find_package(LAPACK REQUIRED)
find_package(Threads REQUIRED)


#------------------------------------------------------------
# Sources
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

set(MUMPS_COMMON_SRCS src/mumps_part9.F
                      src/mumps_common.c
                      src/mumps_ooc_common.F
                      src/mumps_orderings.c
                      src/mumps_size.c
                      src/mumps_io.c
                      src/mumps_io_basic.c
                      src/mumps_io_thread.c
                      src/mumps_io_err.c
                      src/mumps_static_mapping.F
                      src/tools_common_mod.F
                      src/mumps_sol_es.F)

set(MUMPS_C_SRCS src/cmumps_part1.F
                 src/cmumps_part2.F
                 src/cmumps_part3.F
                 src/cmumps_part4.F
                 src/cmumps_part5.F
                 src/cmumps_part6.F
                 src/cmumps_part7.F
                 src/cmumps_part8.F
                 src/cmumps_comm_buffer.F
                 src/cmumps_load.F
                 src/mumps_c.c
                 src/cmumps_ooc_buffer.F
                 src/cmumps_ooc.F
                 src/cmumps_struc_def.F)

set(MUMPS_D_SRCS src/dmumps_part1.F
                 src/dmumps_part2.F
                 src/dmumps_part3.F
                 src/dmumps_part4.F
                 src/dmumps_part5.F
                 src/dmumps_part6.F
                 src/dmumps_part7.F
                 src/dmumps_part8.F
                 src/dmumps_comm_buffer.F
                 src/dmumps_load.F
                 src/mumps_c.c
                 src/dmumps_ooc_buffer.F
                 src/dmumps_ooc.F
                 src/dmumps_struc_def.F)

set(MUMPS_S_SRCS src/smumps_part1.F
                 src/smumps_part2.F
                 src/smumps_part3.F
                 src/smumps_part4.F
                 src/smumps_part5.F
                 src/smumps_part6.F
                 src/smumps_part7.F
                 src/smumps_part8.F
                 src/smumps_comm_buffer.F
                 src/smumps_load.F
                 src/mumps_c.c
                 src/smumps_ooc_buffer.F
                 src/smumps_ooc.F
                 src/smumps_struc_def.F)

set(MUMPS_Z_SRCS src/zmumps_part1.F
                 src/zmumps_part2.F
                 src/zmumps_part3.F
                 src/zmumps_part4.F
                 src/zmumps_part5.F
                 src/zmumps_part6.F
                 src/zmumps_part7.F
                 src/zmumps_part8.F
                 src/zmumps_comm_buffer.F
                 src/zmumps_load.F
                 src/mumps_c.c
                 src/zmumps_ooc_buffer.F
                 src/zmumps_ooc.F
                 src/zmumps_struc_def.F)



# Warnings are disabled
if (MSVC)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /fpp /nologo /reentrancy /fixed /warn:noalignments /Qsave /Qzero /libs:static /threads /traceback /D_CRT_SECURE_NO_WARNINGS /DALLOW_NON_INIT /Dintel_ ")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /nologo /D_CRT_SECURE_NO_WARNINGS /DAdd_ ")
else ()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -w -fcray-pointer -fall-intrinsics -finit-local-zero -DALLOW_NON_INIT -Dintel_ ")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w -DAdd_ ")
endif ()


if(MUMPS_BUILD_SHARED_LIBS)
  add_library(coinmumps SHARED ${MUMPS_D_SRCS} ${MUMPS_COMMON_SRCS} ${MUMPS_LIBSEQ_SRCS})
else()
  add_library(coinmumps STATIC ${MUMPS_D_SRCS} ${MUMPS_COMMON_SRCS} ${MUMPS_LIBSEQ_SRCS})
endif()

if (MSVC)
  target_compile_definitions(coinmumps PRIVATE "/DMUMPS_ARITH=MUMPS_ARITH_d")
else ()
  target_compile_definitions(coinmumps PRIVATE "-DMUMPS_ARITH=MUMPS_ARITH_d")
endif ()
target_link_libraries(coinmumps PUBLIC ${LAPACK_LIBRARIES})

if(THREADS_HAVE_PTHREAD_ARG)
  target_compile_options(coinmumps PUBLIC pthread)
endif()
if(CMAKE_THREAD_LIBS_INIT)
  target_link_libraries(coinmumps PUBLIC ${CMAKE_THREAD_LIBS_INIT})
endif()

target_include_directories(coinmumps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(coinmumps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libseq)


# if (MUMPS_USE_METIS)
#   target_link_libraries(dmumps PUBLIC metis m)
#   target_compile_definitions(dmumps PRIVATE -Dmetis)
# endif ()

## If you want to check dmumps lib works fine.
## See examples/README for instructions running the tests
# add_executable(dsimple_test examples/dsimpletest.F)
# target_link_libraries(dsimple_test dmumps ${LINK_LIBS})



install(TARGETS coinmumps)


# ----------------------------------------------------------------------


# # If you want to build the other mumps libs, e.g. single precision (smumps) or complex mumps (cmumps)
# # copy what is done for dmumps above and adjust it. It is also a good idea to build mumps_common,
# # libseq, and pord as separate libraries in that case. Right now since we only want dmumps we just
# # build it using all need sources to avoid creating unnecessary libraries that are linked just once.

# add_library(seq STATIC ${MUMPS_LIBSEQ_SRCS})

# add_library(pord STATIC ${MUMPS_PORD_SRCS})
# target_include_directories(pord BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/PORD/include)

# add_library(mumps_common STATIC ${MUMPS_COMMON_SRCS})



# add_library(cmumps STATIC ${MUMPS_C_SRCS})
# target_link_libraries(cmumps mumps_common)

# add_library(smumps STATIC ${MUMPS_S_SRCS})
# target_link_libraries(smumps mumps_common)

# add_library(zmumps STATIC ${MUMPS_Z_SRCS})
# target_link_libraries(zmumps mumps_common)

# if (MSVC)
#   target_compile_definitions(smumps PRIVATE "/DMUMPS_ARITH=MUMPS_ARITH_s")
#   target_compile_definitions(dmumps PRIVATE "/DMUMPS_ARITH=MUMPS_ARITH_d")
#   target_compile_definitions(cmumps PRIVATE "/DMUMPS_ARITH=MUMPS_ARITH_c")
#   target_compile_definitions(zmumps PRIVATE "/DMUMPS_ARITH=MUMPS_ARITH_z")
# else ()
#   target_compile_definitions(smumps PRIVATE "-DMUMPS_ARITH=MUMPS_ARITH_s")
#   target_compile_definitions(dmumps PRIVATE "-DMUMPS_ARITH=MUMPS_ARITH_d")
#   target_compile_definitions(cmumps PRIVATE "-DMUMPS_ARITH=MUMPS_ARITH_c")
#   target_compile_definitions(zmumps PRIVATE "-DMUMPS_ARITH=MUMPS_ARITH_z")
# endif ()



