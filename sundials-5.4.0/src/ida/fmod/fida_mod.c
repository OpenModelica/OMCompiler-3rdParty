/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.0
 *
 * This file is not intended to be easily readable and contains a number of
 * coding conventions designed to improve portability and efficiency. Do not make
 * changes to this file unless you know what you are doing--modify the SWIG
 * interface file instead.
 * ----------------------------------------------------------------------------- */

/* ---------------------------------------------------------------
 * Programmer(s): Auto-generated by swig.
 * ---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -------------------------------------------------------------*/

/* -----------------------------------------------------------------------------
 *  This section contains generic SWIG labels for method/variable
 *  declarations/attributes, and other compiler dependent labels.
 * ----------------------------------------------------------------------------- */

/* template workaround for compilers that cannot correctly implement the C++ standard */
#ifndef SWIGTEMPLATEDISAMBIGUATOR
# if defined(__SUNPRO_CC) && (__SUNPRO_CC <= 0x560)
#  define SWIGTEMPLATEDISAMBIGUATOR template
# elif defined(__HP_aCC)
/* Needed even with `aCC -AA' when `aCC -V' reports HP ANSI C++ B3910B A.03.55 */
/* If we find a maximum version that requires this, the test would be __HP_aCC <= 35500 for A.03.55 */
#  define SWIGTEMPLATEDISAMBIGUATOR template
# else
#  define SWIGTEMPLATEDISAMBIGUATOR
# endif
#endif

/* inline attribute */
#ifndef SWIGINLINE
# if defined(__cplusplus) || (defined(__GNUC__) && !defined(__STRICT_ANSI__))
#   define SWIGINLINE inline
# else
#   define SWIGINLINE
# endif
#endif

/* attribute recognised by some compilers to avoid 'unused' warnings */
#ifndef SWIGUNUSED
# if defined(__GNUC__)
#   if !(defined(__cplusplus)) || (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#     define SWIGUNUSED __attribute__ ((__unused__))
#   else
#     define SWIGUNUSED
#   endif
# elif defined(__ICC)
#   define SWIGUNUSED __attribute__ ((__unused__))
# else
#   define SWIGUNUSED
# endif
#endif

#ifndef SWIG_MSC_UNSUPPRESS_4505
# if defined(_MSC_VER)
#   pragma warning(disable : 4505) /* unreferenced local function has been removed */
# endif
#endif

#ifndef SWIGUNUSEDPARM
# ifdef __cplusplus
#   define SWIGUNUSEDPARM(p)
# else
#   define SWIGUNUSEDPARM(p) p SWIGUNUSED
# endif
#endif

/* internal SWIG method */
#ifndef SWIGINTERN
# define SWIGINTERN static SWIGUNUSED
#endif

/* internal inline SWIG method */
#ifndef SWIGINTERNINLINE
# define SWIGINTERNINLINE SWIGINTERN SWIGINLINE
#endif

/* qualifier for exported *const* global data variables*/
#ifndef SWIGEXTERN
# ifdef __cplusplus
#   define SWIGEXTERN extern
# else
#   define SWIGEXTERN
# endif
#endif

/* exporting methods */
#if defined(__GNUC__)
#  if (__GNUC__ >= 4) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#    ifndef GCC_HASCLASSVISIBILITY
#      define GCC_HASCLASSVISIBILITY
#    endif
#  endif
#endif

#ifndef SWIGEXPORT
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   if defined(STATIC_LINKED)
#     define SWIGEXPORT
#   else
#     define SWIGEXPORT __declspec(dllexport)
#   endif
# else
#   if defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#     define SWIGEXPORT __attribute__ ((visibility("default")))
#   else
#     define SWIGEXPORT
#   endif
# endif
#endif

/* calling conventions for Windows */
#ifndef SWIGSTDCALL
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   define SWIGSTDCALL __stdcall
# else
#   define SWIGSTDCALL
# endif
#endif

/* Deal with Microsoft's attempt at deprecating C standard runtime functions */
#if !defined(SWIG_NO_CRT_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_CRT_SECURE_NO_DEPRECATE)
# define _CRT_SECURE_NO_DEPRECATE
#endif

/* Deal with Microsoft's attempt at deprecating methods in the standard C++ library */
#if !defined(SWIG_NO_SCL_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_SCL_SECURE_NO_DEPRECATE)
# define _SCL_SECURE_NO_DEPRECATE
#endif

/* Deal with Apple's deprecated 'AssertMacros.h' from Carbon-framework */
#if defined(__APPLE__) && !defined(__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES)
# define __ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES 0
#endif

/* Intel's compiler complains if a variable which was never initialised is
 * cast to void, which is a common idiom which we use to indicate that we
 * are aware a variable isn't used.  So we just silence that warning.
 * See: https://github.com/swig/swig/issues/192 for more discussion.
 */
#ifdef __INTEL_COMPILER
# pragma warning disable 592
#endif

/*  Errors in SWIG */
#define  SWIG_UnknownError    	   -1
#define  SWIG_IOError        	   -2
#define  SWIG_RuntimeError   	   -3
#define  SWIG_IndexError     	   -4
#define  SWIG_TypeError      	   -5
#define  SWIG_DivisionByZero 	   -6
#define  SWIG_OverflowError  	   -7
#define  SWIG_SyntaxError    	   -8
#define  SWIG_ValueError     	   -9
#define  SWIG_SystemError    	   -10
#define  SWIG_AttributeError 	   -11
#define  SWIG_MemoryError    	   -12
#define  SWIG_NullReferenceError   -13




#include <assert.h>
#define SWIG_exception_impl(DECL, CODE, MSG, RETURNNULL) \
 { printf("In " DECL ": " MSG); assert(0); RETURNNULL; }


#include <stdio.h>
#if defined(_MSC_VER) || defined(__BORLANDC__) || defined(_WATCOM)
# ifndef snprintf
#  define snprintf _snprintf
# endif
#endif


/* Support for the `contract` feature.
 *
 * Note that RETURNNULL is first because it's inserted via a 'Replaceall' in
 * the fortran.cxx file.
 */
#define SWIG_contract_assert(RETURNNULL, EXPR, MSG) \
 if (!(EXPR)) { SWIG_exception_impl("$decl", SWIG_ValueError, MSG, RETURNNULL); } 


#define SWIGVERSION 0x040000 
#define SWIG_VERSION SWIGVERSION


#define SWIG_as_voidptr(a) (void *)((const void *)(a)) 
#define SWIG_as_voidptrptr(a) ((void)SWIG_as_voidptr(*a),(void**)(a)) 


#include "ida/ida.h"
#include "ida/ida_bbdpre.h"
#include "ida/ida_ls.h"


#include <stdlib.h>
#ifdef _MSC_VER
# ifndef strtoull
#  define strtoull _strtoui64
# endif
# ifndef strtoll
#  define strtoll _strtoi64
# endif
#endif


typedef struct {
    void* data;
    size_t size;
} SwigArrayWrapper;


SWIGINTERN SwigArrayWrapper SwigArrayWrapper_uninitialized() {
  SwigArrayWrapper result;
  result.data = NULL;
  result.size = 0;
  return result;
}


#include <string.h>

SWIGEXPORT void * _wrap_FIDACreate() {
  void * fresult ;
  void *result = 0 ;
  
  result = (void *)IDACreate();
  fresult = result;
  return fresult;
}


SWIGEXPORT int _wrap_FIDAInit(void *farg1, IDAResFn farg2, double const *farg3, N_Vector farg4, N_Vector farg5) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  IDAResFn arg2 = (IDAResFn) 0 ;
  realtype arg3 ;
  N_Vector arg4 = (N_Vector) 0 ;
  N_Vector arg5 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (IDAResFn)(farg2);
  arg3 = (realtype)(*farg3);
  arg4 = (N_Vector)(farg4);
  arg5 = (N_Vector)(farg5);
  result = (int)IDAInit(arg1,arg2,arg3,arg4,arg5);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAReInit(void *farg1, double const *farg2, N_Vector farg3, N_Vector farg4) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  N_Vector arg3 = (N_Vector) 0 ;
  N_Vector arg4 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  arg3 = (N_Vector)(farg3);
  arg4 = (N_Vector)(farg4);
  result = (int)IDAReInit(arg1,arg2,arg3,arg4);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASStolerances(void *farg1, double const *farg2, double const *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  realtype arg3 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  arg3 = (realtype)(*farg3);
  result = (int)IDASStolerances(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASVtolerances(void *farg1, double const *farg2, N_Vector farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  N_Vector arg3 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  arg3 = (N_Vector)(farg3);
  result = (int)IDASVtolerances(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAWFtolerances(void *farg1, IDAEwtFn farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  IDAEwtFn arg2 = (IDAEwtFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (IDAEwtFn)(farg2);
  result = (int)IDAWFtolerances(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDACalcIC(void *farg1, int const *farg2, double const *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  realtype arg3 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  arg3 = (realtype)(*farg3);
  result = (int)IDACalcIC(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetNonlinConvCoefIC(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)IDASetNonlinConvCoefIC(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetMaxNumStepsIC(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetMaxNumStepsIC(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetMaxNumJacsIC(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetMaxNumJacsIC(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetMaxNumItersIC(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetMaxNumItersIC(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetLineSearchOffIC(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetLineSearchOffIC(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetStepToleranceIC(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)IDASetStepToleranceIC(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetMaxBacksIC(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetMaxBacksIC(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetErrHandlerFn(void *farg1, IDAErrHandlerFn farg2, void *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  IDAErrHandlerFn arg2 = (IDAErrHandlerFn) 0 ;
  void *arg3 = (void *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (IDAErrHandlerFn)(farg2);
  arg3 = (void *)(farg3);
  result = (int)IDASetErrHandlerFn(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetErrFile(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  FILE *arg2 = (FILE *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (FILE *)(farg2);
  result = (int)IDASetErrFile(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetUserData(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  void *arg2 = (void *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (void *)(farg2);
  result = (int)IDASetUserData(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetMaxOrd(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetMaxOrd(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetMaxNumSteps(void *farg1, long const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long)(*farg2);
  result = (int)IDASetMaxNumSteps(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetInitStep(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)IDASetInitStep(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetMaxStep(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)IDASetMaxStep(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetStopTime(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)IDASetStopTime(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetNonlinConvCoef(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)IDASetNonlinConvCoef(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetMaxErrTestFails(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetMaxErrTestFails(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetMaxNonlinIters(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetMaxNonlinIters(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetMaxConvFails(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetMaxConvFails(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetSuppressAlg(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetSuppressAlg(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetId(void *farg1, N_Vector farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector arg2 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector)(farg2);
  result = (int)IDASetId(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetConstraints(void *farg1, N_Vector farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector arg2 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector)(farg2);
  result = (int)IDASetConstraints(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetNonlinearSolver(void *farg1, SUNNonlinearSolver farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  SUNNonlinearSolver arg2 = (SUNNonlinearSolver) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (SUNNonlinearSolver)(farg2);
  result = (int)IDASetNonlinearSolver(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDARootInit(void *farg1, int const *farg2, IDARootFn farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  IDARootFn arg3 = (IDARootFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  arg3 = (IDARootFn)(farg3);
  result = (int)IDARootInit(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetRootDirection(void *farg1, int *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int *arg2 = (int *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int *)(farg2);
  result = (int)IDASetRootDirection(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetNoInactiveRootWarn(void *farg1) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  result = (int)IDASetNoInactiveRootWarn(arg1);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASolve(void *farg1, double const *farg2, double *farg3, N_Vector farg4, N_Vector farg5, int const *farg6) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  realtype *arg3 = (realtype *) 0 ;
  N_Vector arg4 = (N_Vector) 0 ;
  N_Vector arg5 = (N_Vector) 0 ;
  int arg6 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  arg3 = (realtype *)(farg3);
  arg4 = (N_Vector)(farg4);
  arg5 = (N_Vector)(farg5);
  arg6 = (int)(*farg6);
  result = (int)IDASolve(arg1,arg2,arg3,arg4,arg5,arg6);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAComputeY(void *farg1, N_Vector farg2, N_Vector farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector arg2 = (N_Vector) 0 ;
  N_Vector arg3 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector)(farg2);
  arg3 = (N_Vector)(farg3);
  result = (int)IDAComputeY(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAComputeYp(void *farg1, N_Vector farg2, N_Vector farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector arg2 = (N_Vector) 0 ;
  N_Vector arg3 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector)(farg2);
  arg3 = (N_Vector)(farg3);
  result = (int)IDAComputeYp(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetDky(void *farg1, double const *farg2, int const *farg3, N_Vector farg4) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int arg3 ;
  N_Vector arg4 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  arg3 = (int)(*farg3);
  arg4 = (N_Vector)(farg4);
  result = (int)IDAGetDky(arg1,arg2,arg3,arg4);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetWorkSpace(void *farg1, long *farg2, long *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  long *arg3 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  arg3 = (long *)(farg3);
  result = (int)IDAGetWorkSpace(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumSteps(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumSteps(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumResEvals(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumResEvals(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumLinSolvSetups(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumLinSolvSetups(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumErrTestFails(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumErrTestFails(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumBacktrackOps(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumBacktrackOps(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetConsistentIC(void *farg1, N_Vector farg2, N_Vector farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector arg2 = (N_Vector) 0 ;
  N_Vector arg3 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector)(farg2);
  arg3 = (N_Vector)(farg3);
  result = (int)IDAGetConsistentIC(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetLastOrder(void *farg1, int *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int *arg2 = (int *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int *)(farg2);
  result = (int)IDAGetLastOrder(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetCurrentOrder(void *farg1, int *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int *arg2 = (int *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int *)(farg2);
  result = (int)IDAGetCurrentOrder(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetCurrentCj(void *farg1, double *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype *arg2 = (realtype *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype *)(farg2);
  result = (int)IDAGetCurrentCj(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetCurrentY(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector *arg2 = (N_Vector *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector *)(farg2);
  result = (int)IDAGetCurrentY(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetCurrentYp(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector *arg2 = (N_Vector *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector *)(farg2);
  result = (int)IDAGetCurrentYp(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetActualInitStep(void *farg1, double *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype *arg2 = (realtype *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype *)(farg2);
  result = (int)IDAGetActualInitStep(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetLastStep(void *farg1, double *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype *arg2 = (realtype *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype *)(farg2);
  result = (int)IDAGetLastStep(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetCurrentStep(void *farg1, double *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype *arg2 = (realtype *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype *)(farg2);
  result = (int)IDAGetCurrentStep(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetCurrentTime(void *farg1, double *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype *arg2 = (realtype *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype *)(farg2);
  result = (int)IDAGetCurrentTime(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetTolScaleFactor(void *farg1, double *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype *arg2 = (realtype *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype *)(farg2);
  result = (int)IDAGetTolScaleFactor(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetErrWeights(void *farg1, N_Vector farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector arg2 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector)(farg2);
  result = (int)IDAGetErrWeights(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetEstLocalErrors(void *farg1, N_Vector farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector arg2 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector)(farg2);
  result = (int)IDAGetEstLocalErrors(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumGEvals(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumGEvals(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetRootInfo(void *farg1, int *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int *arg2 = (int *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int *)(farg2);
  result = (int)IDAGetRootInfo(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetIntegratorStats(void *farg1, long *farg2, long *farg3, long *farg4, long *farg5, int *farg6, int *farg7, double *farg8, double *farg9, double *farg10, double *farg11) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  long *arg3 = (long *) 0 ;
  long *arg4 = (long *) 0 ;
  long *arg5 = (long *) 0 ;
  int *arg6 = (int *) 0 ;
  int *arg7 = (int *) 0 ;
  realtype *arg8 = (realtype *) 0 ;
  realtype *arg9 = (realtype *) 0 ;
  realtype *arg10 = (realtype *) 0 ;
  realtype *arg11 = (realtype *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  arg3 = (long *)(farg3);
  arg4 = (long *)(farg4);
  arg5 = (long *)(farg5);
  arg6 = (int *)(farg6);
  arg7 = (int *)(farg7);
  arg8 = (realtype *)(farg8);
  arg9 = (realtype *)(farg9);
  arg10 = (realtype *)(farg10);
  arg11 = (realtype *)(farg11);
  result = (int)IDAGetIntegratorStats(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNonlinearSystemData(void *farg1, double *farg2, void *farg3, void *farg4, void *farg5, void *farg6, void *farg7, double *farg8, void *farg9) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype *arg2 = (realtype *) 0 ;
  N_Vector *arg3 = (N_Vector *) 0 ;
  N_Vector *arg4 = (N_Vector *) 0 ;
  N_Vector *arg5 = (N_Vector *) 0 ;
  N_Vector *arg6 = (N_Vector *) 0 ;
  N_Vector *arg7 = (N_Vector *) 0 ;
  realtype *arg8 = (realtype *) 0 ;
  void **arg9 = (void **) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype *)(farg2);
  arg3 = (N_Vector *)(farg3);
  arg4 = (N_Vector *)(farg4);
  arg5 = (N_Vector *)(farg5);
  arg6 = (N_Vector *)(farg6);
  arg7 = (N_Vector *)(farg7);
  arg8 = (realtype *)(farg8);
  arg9 = (void **)(farg9);
  result = (int)IDAGetNonlinearSystemData(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumNonlinSolvIters(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumNonlinSolvIters(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumNonlinSolvConvFails(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumNonlinSolvConvFails(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNonlinSolvStats(void *farg1, long *farg2, long *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  long *arg3 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  arg3 = (long *)(farg3);
  result = (int)IDAGetNonlinSolvStats(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT SwigArrayWrapper _wrap_FIDAGetReturnFlagName(long const *farg1) {
  SwigArrayWrapper fresult ;
  long arg1 ;
  char *result = 0 ;
  
  arg1 = (long)(*farg1);
  result = (char *)IDAGetReturnFlagName(arg1);
  fresult.size = strlen((const char*)(result));
  fresult.data = (char *)(result);
  return fresult;
}


SWIGEXPORT void _wrap_FIDAFree(void *farg1) {
  void **arg1 = (void **) 0 ;
  
  arg1 = (void **)(farg1);
  IDAFree(arg1);
}


SWIGEXPORT int _wrap_FIDASetJacTimesResFn(void *farg1, IDAResFn farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  IDAResFn arg2 = (IDAResFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (IDAResFn)(farg2);
  result = (int)IDASetJacTimesResFn(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDABBDPrecInit(void *farg1, int64_t const *farg2, int64_t const *farg3, int64_t const *farg4, int64_t const *farg5, int64_t const *farg6, double const *farg7, IDABBDLocalFn farg8, IDABBDCommFn farg9) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  sunindextype arg2 ;
  sunindextype arg3 ;
  sunindextype arg4 ;
  sunindextype arg5 ;
  sunindextype arg6 ;
  realtype arg7 ;
  IDABBDLocalFn arg8 = (IDABBDLocalFn) 0 ;
  IDABBDCommFn arg9 = (IDABBDCommFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (sunindextype)(*farg2);
  arg3 = (sunindextype)(*farg3);
  arg4 = (sunindextype)(*farg4);
  arg5 = (sunindextype)(*farg5);
  arg6 = (sunindextype)(*farg6);
  arg7 = (realtype)(*farg7);
  arg8 = (IDABBDLocalFn)(farg8);
  arg9 = (IDABBDCommFn)(farg9);
  result = (int)IDABBDPrecInit(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDABBDPrecReInit(void *farg1, int64_t const *farg2, int64_t const *farg3, double const *farg4) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  sunindextype arg2 ;
  sunindextype arg3 ;
  realtype arg4 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (sunindextype)(*farg2);
  arg3 = (sunindextype)(*farg3);
  arg4 = (realtype)(*farg4);
  result = (int)IDABBDPrecReInit(arg1,arg2,arg3,arg4);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDABBDPrecGetWorkSpace(void *farg1, long *farg2, long *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  long *arg3 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  arg3 = (long *)(farg3);
  result = (int)IDABBDPrecGetWorkSpace(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDABBDPrecGetNumGfnEvals(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDABBDPrecGetNumGfnEvals(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetLinearSolver(void *farg1, SUNLinearSolver farg2, SUNMatrix farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  SUNLinearSolver arg2 = (SUNLinearSolver) 0 ;
  SUNMatrix arg3 = (SUNMatrix) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (SUNLinearSolver)(farg2);
  arg3 = (SUNMatrix)(farg3);
  result = (int)IDASetLinearSolver(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetJacFn(void *farg1, IDALsJacFn farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  IDALsJacFn arg2 = (IDALsJacFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (IDALsJacFn)(farg2);
  result = (int)IDASetJacFn(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetPreconditioner(void *farg1, IDALsPrecSetupFn farg2, IDALsPrecSolveFn farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  IDALsPrecSetupFn arg2 = (IDALsPrecSetupFn) 0 ;
  IDALsPrecSolveFn arg3 = (IDALsPrecSolveFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (IDALsPrecSetupFn)(farg2);
  arg3 = (IDALsPrecSolveFn)(farg3);
  result = (int)IDASetPreconditioner(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetJacTimes(void *farg1, IDALsJacTimesSetupFn farg2, IDALsJacTimesVecFn farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  IDALsJacTimesSetupFn arg2 = (IDALsJacTimesSetupFn) 0 ;
  IDALsJacTimesVecFn arg3 = (IDALsJacTimesVecFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (IDALsJacTimesSetupFn)(farg2);
  arg3 = (IDALsJacTimesVecFn)(farg3);
  result = (int)IDASetJacTimes(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetEpsLin(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)IDASetEpsLin(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetLSNormFactor(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)IDASetLSNormFactor(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetLinearSolutionScaling(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)IDASetLinearSolutionScaling(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDASetIncrementFactor(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)IDASetIncrementFactor(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetLinWorkSpace(void *farg1, long *farg2, long *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  long *arg3 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  arg3 = (long *)(farg3);
  result = (int)IDAGetLinWorkSpace(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumJacEvals(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumJacEvals(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumPrecEvals(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumPrecEvals(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumPrecSolves(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumPrecSolves(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumLinIters(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumLinIters(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumLinConvFails(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumLinConvFails(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumJTSetupEvals(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumJTSetupEvals(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumJtimesEvals(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumJtimesEvals(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetNumLinResEvals(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetNumLinResEvals(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FIDAGetLastLinFlag(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)IDAGetLastLinFlag(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT SwigArrayWrapper _wrap_FIDAGetLinReturnFlagName(long const *farg1) {
  SwigArrayWrapper fresult ;
  long arg1 ;
  char *result = 0 ;
  
  arg1 = (long)(*farg1);
  result = (char *)IDAGetLinReturnFlagName(arg1);
  fresult.size = strlen((const char*)(result));
  fresult.data = (char *)(result);
  return fresult;
}


