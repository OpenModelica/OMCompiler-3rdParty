cmake_minimum_required(VERSION 3.14)
project(Lis)

# If we are compiling for an x86 CPU set HAS_X87_FPU.
if(OMC_TARGET_ARCH MATCHES "x86")
    set (HAS_X87_FPU TRUE)
else()
    set (HAS_X87_FPU FALSE)
endif()

check_type_size(double SIZEOF_DOUBLE)
check_type_size(float SIZEOF_FLOAT)
check_type_size(int SIZEOF_INT)
check_type_size(long SIZEOF_LONG)
check_type_size("long double" SIZEOF_LONG_DOUBLE)
check_type_size("long long" SIZEOF_LONG_LONG)
check_type_size(size_t SIZEOF_SIZE_T)
set(SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})


# Check if we have the headers needed
set(HEADERS_TO_CHECK dlfcn.h inttypes.h malloc.h memory.h stdint.h stdlib.h strings.h
                     string.h sys/stat.h sys/time.h sys/types.h unistd.h sys/mount.h)
# The quote is needed here to evaluate
omc_check_headers_exist_and_define_each("${HEADERS_TO_CHECK}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/lis_config.h.in.cmake ${CMAKE_CURRENT_SOURCE_DIR}/include/lis_config.h)

add_subdirectory(src)
