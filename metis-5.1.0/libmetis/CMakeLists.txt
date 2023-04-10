
# Add this directory for internal users.
include_directories(.)
# Find sources.
file(GLOB metis_sources *.c)
# Build libmetis.
add_library(metis ${METIS_LIBRARY_TYPE} ${GKlib_sources} ${metis_sources})

if(UNIX)
  target_link_libraries(metis m)
elseif(MSVC)
  target_link_libraries(metis PUBLIC omc::3rd::regex)
endif()

if(METIS_INSTALL)
  install(TARGETS metis DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
