find_package (PkgConfig QUIET)
if (PkgConfig_FOUND)
  pkg_check_modules (PC_LIBFFI QUIET libffi)
endif ()

# Find the headers and library
find_path(
  LIBFFI_INCLUDE_DIR
  NAMES "ffi.h"
  HINTS "${PC_LIBFFI_INCLUDE_DIRS}")

find_library(
  LIBFFI_LIBRARY
  NAMES "ffi"
  HINTS "${PC_LIBFFI_LIBRARY_DIRS}")

set (LIBFFI_LIBRARIES ${LIBFFI_LIBRARY})
set (LIBFFI_INCLUDE_DIRS ${LIBFFI_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibFFI DEFAULT_MSG LIBFFI_LIBRARY LIBFFI_INCLUDE_DIRS)

mark_as_advanced (
  LIBFFI_LIBRARY
  LIBFFI_LIBRARIES
  LIBFFI_INCLUDE_DIR
  LIBFFI_INCLUDE_DIRS)
