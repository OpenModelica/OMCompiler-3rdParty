# The name of our project is "CMINPACK". CMakeLists files in this project can
# refer to the root source directory of the project as ${CMINPACK_SOURCE_DIR} and
# to the root binary directory of the project as ${CMINPACK_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (CMINPACK)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)

# Use whatever is set in CMAKE_INSTALL_LIBDIR instead of computing a directory
# internally. For OMC this is set from outside.
set(CMINPACK_LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})

include(${PROJECT_SOURCE_DIR}/cmake/cminpack_utils.cmake)
# Set version and OS-specific settings
set(CMINPACK_VERSION 1.3.4 CACHE STRING "CMinpack version")
set(CMINPACK_SOVERSION 1 CACHE STRING "CMinpack API version")
DISSECT_VERSION()
GET_OS_INFO()

# Add an "uninstall" target
CONFIGURE_FILE ("${PROJECT_SOURCE_DIR}/cmake/uninstall_target.cmake.in"
    "${PROJECT_BINARY_DIR}/uninstall_target.cmake" IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET (uninstall "${CMAKE_COMMAND}" -P
    "${PROJECT_BINARY_DIR}/uninstall_target.cmake")

enable_testing()

if (OS_LINUX)
  option (CMINPACK_USE_FPIC "Use the -fPIC compiler flag." ON)
else (OS_LINUX)
  option (CMINPACK_USE_FPIC "Use the -fPIC compiler flag." OFF)
endif (OS_LINUX)

option (CMINPACK_BUILD_SHARED_LIBS "Build shared libraries instead of static." OFF)
if (CMINPACK_BUILD_SHARED_LIBS)
  set (CMINPACK_LIBRARY_BUILD_TYPE SHARED)
  message (STATUS "Building shared libraries.")
else ()
  message (STATUS "Building static libraries.")
  set (CMINPACK_LIBRARY_BUILD_TYPE STATIC)
  if(WIN32)
    add_definitions(-DCMINPACK_NO_DLL)
  endif(WIN32)
endif ()

#set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../build)

add_subdirectory (cmake)
#add_subdirectory (examples)

set (cminpack_srcs
  cminpack.h cminpackP.h
  chkder.c  enorm.c   hybrd1.c  hybrj.c   lmdif1.c  lmstr1.c  qrfac.c   r1updt.c
  dogleg.c  fdjac1.c  hybrd.c   lmder1.c  lmdif.c   lmstr.c   qrsolv.c  rwupdt.c
  dpmpar.c  fdjac2.c  hybrj1.c  lmder.c   lmpar.c   qform.c   r1mpyq.c  covar.c covar1.c
  minpack.h
  chkder_.c enorm_.c  hybrd1_.c hybrj_.c  lmdif1_.c lmstr1_.c qrfac_.c  r1updt_.c
  dogleg_.c fdjac1_.c hybrd_.c  lmder1_.c lmdif_.c  lmstr_.c  qrsolv_.c rwupdt_.c
  dpmpar_.c fdjac2_.c hybrj1_.c lmder_.c  lmpar_.c  qform_.c  r1mpyq_.c covar_.c
  )
set (cminpack_hdrs
    cminpack.h minpack.h)

add_library (cminpack ${CMINPACK_LIBRARY_BUILD_TYPE} ${cminpack_srcs})

install (TARGETS cminpack
   LIBRARY DESTINATION ${CMINPACK_LIB_INSTALL_DIR}
   ARCHIVE DESTINATION ${CMINPACK_LIB_INSTALL_DIR}
   RUNTIME DESTINATION bin)
install (FILES ${cminpack_hdrs} DESTINATION ${CMINPACK_INCLUDE_INSTALL_DIR})

if (CMINPACK_USE_FPIC AND NOT CMINPACK_BUILD_SHARED_LIBS)
  set_target_properties (cminpack PROPERTIES COMPILE_FLAGS -fPIC)
endif ()

set_target_properties(cminpack PROPERTIES VERSION ${CMINPACK_VERSION} SOVERSION ${CMINPACK_SOVERSION})

