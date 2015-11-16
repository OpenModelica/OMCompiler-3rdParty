project(umfpack_package)
cmake_minimum_required(VERSION 2.8.9)
include(ExternalProject)
include(Helpers.cmake)
set(SS_VERSION "3.0.0")

set(LIB_OMC "lib/omc" CACHE STRING "library directory")

#for shared, at least in cmake 2.8.12:
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(DOWNLOAD_SUITESPARSE "download suitesparse package, version ${SS_VERSION}" OFF)
option(WITH_BLAS "enable blas support" OFF)
option(WITH_CHOLMOD "enable cholmod support" OFF)

if(DOWNLOAD_SUITESPARSE)
  if(NOT EXISTS "SuiteSparse.tar.gz")
    file(DOWNLOAD "http://www.cise.ufl.edu/research/sparse/SuiteSparse/SuiteSparse-${SS_VERSION}.tar.gz"
    ./SuiteSparse.tar.gz SHOW_PROGRESS)
  endif()
  execute_process(COMMAND cmake -E tar xzf "SuiteSparse.tar.gz")
endif()

find_file(SS_INSTALL README.txt DOC "README.txt in SuiteSparse folder"
  HINT "${CMAKE_BINARY_DIR}/.." NO_DEFAULT_PATH)
if(SS_INSTALL)
  get_filename_component(SS_DIR ${SS_INSTALL} PATH)
  if(NOT EXISTS "${SS_DIR}/README.txt")
    message(FATAL_ERROR "please give absolute path to SuiteSparse_install.m")
  endif()
else()
  set(SS_DIR "SuiteSparse-NOTFOUND")
endif()

if(SS_DIR)
  if(NOT WITH_BLAS)
    add_definitions(-DNBLAS)
  endif()

  if(EXISTS "${SS_DIR}/SuiteSparse_config")
    include_directories(${SS_DIR}/SuiteSparse_config)
    add_library(SuiteSparse_config ${SS_DIR}/SuiteSparse_config/SuiteSparse_config.c)
    install(FILES ${SS_DIR}/SuiteSparse_config/SuiteSparse_config.h
    DESTINATION include/omc/c/suitesparse/Include)
    install(TARGETS SuiteSparse_config DESTINATION ${LIB_OMC})
  elseif(EXISTS "${SS_DIR}/UFconfig")
    include_directories(${SS_DIR}/UFconfig)
    install(FILES ${SS_DIR}/UFconfig/UFconfig.h
    DESTINATION include/omc/c/suitesparse/Include)
  else()
    message(ERROR "unknown suitesparse configuration")
  endif()

  if(WITH_CHOLMOD)
    include_directories(${SS_DIR}/CHOLMOD/Include)
    set(COLAMD_DIR ${SS_DIR}/COLAMD)
    include_directories(${COLAMD_DIR}/Include)

    #TODO: with CAMD and METIS
    add_definitions(-DNCAMD -DNPARTITION)

    #add colamd
    add_library(colamd_object OBJECT "${COLAMD_DIR}/Source/colamd_global.c"
                "${COLAMD_DIR}/Source/colamd.c")
    add_library(colamd_object_l OBJECT "${COLAMD_DIR}/Source/colamd.c")
    get_target_property(FLAGS colamd_object_l COMPILE_DEFINITIONS)
    list(APPEND FLAGS "DLONG")
    set_target_properties(colamd_object_l PROPERTIES COMPILE_DEFINITIONS ${FLAGS})
    add_library(colamd $<TARGET_OBJECTS:colamd_object> $<TARGET_OBJECTS:colamd_object_l>)

    install(TARGETS colamd DESTINATION ${LIB_OMC})
    install(FILES ${COLAMD_DIR}/Include/colamd.h DESTINATION Include)

    set(CHOL_DIRS Core Check Cholesky MatrixOps Partition Modify Supernodal)
    unset(CHOL_SRC_FILES)
    foreach(curDir IN LISTS CHOL_DIRS)
      file(GLOB curFiles "${SS_DIR}/CHOLMOD/${curDir}/*.c")
      foreach(sFile IN LISTS curFiles)
        #get rid of those template files
        get_filename_component(TEST_FOR_T ${sFile} NAME_WE)
        if(${TEST_FOR_T} MATCHES "t_cholmod")
        else()
           list(APPEND CHOL_SRC_FILES ${sFile})
        endif()
      endforeach()
    endforeach()

    add_library(cholmod_object OBJECT ${CHOL_SRC_FILES})
    add_library(cholmod_object_l OBJECT ${CHOL_SRC_FILES})
    get_target_property(FLAGS cholmod_object_l COMPILE_DEFINITIONS)
    list(APPEND FLAGS "DLONG")
    set_target_properties(cholmod_object_l PROPERTIES COMPILE_DEFINITIONS ${FLAGS})
    add_library(cholmod $<TARGET_OBJECTS:cholmod_object> $<TARGET_OBJECTS:cholmod_object_l>)

    install(TARGETS cholmod DESTINATION ${LIB_OMC})
    #collect all cholmod includes
    file(GLOB cholH "${SS_DIR}/CHOLMOD/Include/cholmod*.h")
    foreach(curF IN LISTS cholH)
      if(${curF} MATCHES "cholmod_internal.h")
      else()
        install(FILES "${curF}" DESTINATION "Include")
      endif()
    endforeach()
  else()
    add_definitions(-DNCHOLMOD)
  endif()

  #collect and build the AMD part
  set(AMD_DIR "${SS_DIR}/AMD")
  include_directories(${AMD_DIR}/Include)
  include_directories(${AMD_DIR}/Source)
  install(DIRECTORY ${AMD_DIR}/Include
    DESTINATION include/omc/c/suitesparse)
  file(GLOB curFiles "${AMD_DIR}/Source/*.c")
  foreach(sFile IN LISTS curFiles)
    if(${sFile} MATCHES "amd_global.c")
      add_library(amd_global OBJECT ${sFile})
      set(AMD_OBJECTS "\$<TARGET_OBJECTS:amd_global>")
    else()
      list(APPEND AMD_SRC_FILES ${sFile})
    endif()
  endforeach()
  set(modes "DINT" "DLONG")
  foreach(mode IN LISTS modes)
    add_library(amd_${mode} OBJECT ${AMD_SRC_FILES})
    get_target_property(FLAGS amd_${mode} COMPILE_DEFINITIONS)
    list(APPEND FLAGS ${mode})
    set_target_properties(amd_${mode} PROPERTIES COMPILE_DEFINITIONS "${FLAGS}")
    set(AMD_OBJECTS "${AMD_OBJECTS};\$<TARGET_OBJECTS:amd_${mode}>")
  endforeach()
  add_library(amd ${AMD_OBJECTS})
  install(TARGETS amd DESTINATION ${LIB_OMC})

  #collect and build the BTF part
  set(BTF_DIR "${SS_DIR}/BTF")
  include_directories(${BTF_DIR}/Include)
  include_directories(${BTF_DIR}/Source)
  #install(DIRECTORY ${BTF_DIR}/Include DESTINATION include/omc/c/suitesparse)
  file(GLOB curFiles "${BTF_DIR}/Source/*.c")
  foreach(sFile IN LISTS curFiles)
    list(APPEND BTF_SRC_FILES ${sFile})
  endforeach()
  set(modes "DINT" "DLONG")
  foreach(mode IN LISTS modes)
    add_library(btf_object_${mode} OBJECT ${BTF_SRC_FILES})
    get_target_property(FLAGS btf_object_${mode} COMPILE_DEFINITIONS)
    list(APPEND FLAGS ${mode})
    set_target_properties(btf_object_${mode} PROPERTIES COMPILE_DEFINITIONS "${FLAGS}")
    set(BTF_OBJECTS "${BTF_OBJECTS};\$<TARGET_OBJECTS:btf_object_${mode}>")
  endforeach()
  add_library(btf ${BTF_OBJECTS})
  install(TARGETS btf DESTINATION ${LIB_OMC})
  install(FILES ${BTF_DIR}/Include/btf.h DESTINATION include/omc/c/suitesparse/Include)

  #add colamd
  set(COLAMD_DIR "${SS_DIR}/COLAMD")
  include_directories(${COLAMD_DIR}/Include)
  include_directories(${COLAMD_DIR}/Source)
  add_library(colamd_object OBJECT "${COLAMD_DIR}/Source/colamd_global.c"
              "${COLAMD_DIR}/Source/colamd.c")
  add_library(colamd_object_l OBJECT "${COLAMD_DIR}/Source/colamd.c")
  set_target_properties(colamd_object PROPERTIES COMPILE_DEFINITIONS "DINT")
  set_target_properties(colamd_object_l PROPERTIES COMPILE_DEFINITIONS "DLONG")
  add_library(colamd $<TARGET_OBJECTS:colamd_object> $<TARGET_OBJECTS:colamd_object_l>)
  install(TARGETS colamd DESTINATION ${LIB_OMC})
  install(FILES ${COLAMD_DIR}/Include/colamd.h DESTINATION include/omc/c/suitesparse/Include)

  #add klu
  set(KLU_DIR "${SS_DIR}/KLU")
  include_directories(${KLU_DIR}/Include)
  include_directories(${KLU_DIR}/Source)
  file(GLOB curFiles "${KLU_DIR}/Source/*.c")
  foreach(sFile IN LISTS curFiles)
          list(APPEND KLU_SRC_FILES ${sFile})
  endforeach()
  set(modes "DINT" "DLONG")
  foreach(mode IN LISTS modes)
    add_library(klu_object_${mode} OBJECT ${KLU_SRC_FILES})
    get_target_property(FLAGS klu_object_${mode} COMPILE_DEFINITIONS)
    list(APPEND FLAGS ${mode})
    set_target_properties(klu_object_${mode} PROPERTIES COMPILE_DEFINITIONS "${FLAGS}")
    set(KLU_OBJECTS "${KLU_OBJECTS};\$<TARGET_OBJECTS:klu_object_${mode}>")
  endforeach()
  add_library(klu ${KLU_OBJECTS})
  install(TARGETS klu DESTINATION ${LIB_OMC})
  install(FILES ${KLU_DIR}/Include/klu.h DESTINATION include/omc/c/suitesparse/Include)

  include_directories(${SS_DIR}/UMFPACK/Include)
  include_directories(${SS_DIR}/UMFPACK/Source)
  install(DIRECTORY ${SS_DIR}/UMFPACK/Include
    DESTINATION include/omc/c/suitesparse)
  #collect the umfpack files
  #umfpack structure:
  #   *global, -: umfpack_(global | timer | tictoc).c
  #   *analyze, (DINT, DLONG): umf_(analyze | apply_order | colamd | cholmo |d free | fsize | is_permutation | malloc | realloc | report_perm | singletons).c
  #   *hsolve, (DINT, DLONG, ZINT, ZLONG) | (- | CONJUGATE_SOLVE) : umf_(ltsolve | utsolve).c
  #   *triplet, (DINT, DLONG, ZINT, ZLONG) | (- | DO_MAP | DO_VALUES) : umf_triplet.c
  #   *assemble, (DINT, DLONG, ZINT, ZLONG) | (-, FIXQ) : umf_assemble.c
  #   *store, (DINT, DLONG, ZINT, ZLONG) | (-, DROP) : umf_store_lu.c
  #   *solve, (DINT, DLONG, ZINT, ZLONG) | (-, WSOLVE) : umfpack_solve.c
  #   *others, (DINT, DLONG, ZINT, ZLONG) : *.c

  set(umfpack_dismiss_files umf_multicompile)
  set(umfpack_global_files umfpack_global umfpack_timer umfpack_tictoc)
  set(umfpack_global_groups "-")

  prepandAll(umfpack_analyze_files "umf_" analyze apply_order colamd cholmod free fsize is_permutation malloc realloc report_perm singletons)
  set(umfpack_analyze_groups DINT DLONG)
#  addToFlagGroup("umfpack_olib" "${umfpack_analyze_files}" "DINT;DLONG")

  set(umfpack_hsolve_files umf_ltsolve umf_utsolve)
  getAllGroups(umfpack_hsolve_groups flags_umfpack_olib "DINT;DLONG;ZINT;ZLONG" "-;CONJUGATE_SOLVE")

  set(umfpack_triplet_files umf_triplet)
  getAllGroups(umfpack_triplet_groups flags_umfpack_olib "DINT;DLONG;ZINT;ZLONG" "-;DO_MAP;DO_VALUES;DO_MAP,DO_VALUES")

  set(umfpack_assemble_files umf_assemble)
  getAllGroups(umfpack_assemble_groups flags_umfpack_olib "DINT;DLONG;ZINT;ZLONG" "-;FIXQ")

  set(umfpack_store_files umf_store_lu)
  getAllGroups(umfpack_store_groups flags_umfpack_olib "DINT;DLONG;ZINT;ZLONG" "-;DROP")

  set(umfpack_solve_files umfpack_solve)
  getAllGroups(umfpack_solve_groups flags_umfpack_olib "DINT;DLONG;ZINT;ZLONG" "-;WSOLVE")

  set(umfpack_others_files "*")
  set(umfpack_others_groups "DINT;DLONG;ZINT;ZLONG")

  processGroups(umfpack "global;analyze;hsolve;triplet;assemble;store;solve" "${SS_DIR}/UMFPACK/Source/")
  unset(UMFPACK_OBJECTS )
  foreach( olib IN LISTS list_umfpack_olib)
    add_library(${olib} OBJECT ${${olib}})
    list(LENGTH ${olib} olib_size)
    #message("files ${olib}: ${olib_size}")
    #message("flags_${olib}: ${flags_${olib}}")
    list(LENGTH flags_${olib} olib_size)
    if(${olib_size} GREATER 0)
      get_target_property(FLAGS ${olib} COMPILE_DEFINITIONS)
      foreach(f IN LISTS flags_${olib})
        list(APPEND FLAGS ${f})
      endforeach()
      set_target_properties(${olib} PROPERTIES COMPILE_DEFINITIONS "${FLAGS}")
    endif()
    set(UMFPACK_OBJECTS "${UMFPACK_OBJECTS};\$<TARGET_OBJECTS:${olib}>")
  endforeach()
  add_library(umfpack ${UMFPACK_OBJECTS})
  install(TARGETS umfpack DESTINATION ${LIB_OMC})


  set(PACKAGE_TITLE "small umfpack")
  set(PACKAGE_DESCRIPTION "reduced and selfcontained umfpack. fullversion is available at http://www.cise.ufl.edu/research/sparse/SuiteSparse/")
  set(CPACK_PACKAGE_NAME "small_umfpack")
  set(CPACK_GENERATOR "TGZ" "TBZ2" "ZIP")
  include(CPack)
else()
  message("could not find the suitesparse sources, so won't build anything")
endif()
