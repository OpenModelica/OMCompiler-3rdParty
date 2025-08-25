/*
 *
 *  This file is part of MUMPS 5.8.1, released
 *  on Wed Jul 30 16:49:18 UTC 2025
 *
 *
 *  Copyright 1991-2025 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
#ifndef MUMPS_SCOTCH64_H
#define MUMPS_SCOTCH64_H
#include "mumps_common.h"
#if defined(scotch) || defined(ptscotch)
#include "scotch.h"
/* Instead of the prototypes below, one could include esmumps.h,
 * when provided in include directory of scotch installation */
#if ((SCOTCH_VERSION == 6) && (SCOTCH_RELEASE >= 1)) || (SCOTCH_VERSION >= 7)
/* esmumpsv prototype with weights of nodes in the graph are used on entry (nv) */
MUMPS_INT esmumpsv( const MUMPS_INT8 n, const MUMPS_INT8 iwlen, MUMPS_INT8 * const pe, const MUMPS_INT8 pfree,
             MUMPS_INT8 * const len, MUMPS_INT8 * const iw, MUMPS_INT8 * const nv, MUMPS_INT8 * const elen,
             MUMPS_INT8 * const last);
#endif
/* esmumps prototype (weights of nodes not used on entry) */
MUMPS_INT esmumps( const MUMPS_INT8 n, const MUMPS_INT8 iwlen, MUMPS_INT8 * const pe, const MUMPS_INT8 pfree,
             MUMPS_INT8 * const len, MUMPS_INT8 * const iw, MUMPS_INT8 * const nv, MUMPS_INT8 * const elen,
             MUMPS_INT8 * const last);
#if ((SCOTCH_VERSION == 7) && (SCOTCH_RELEASE >= 1)) || (SCOTCH_VERSION >= 8)
/* prototypes with contexts */
MUMPS_INT esmumpsvc( const MUMPS_INT8 n, const MUMPS_INT8 iwlen, MUMPS_INT8 * const pe, const MUMPS_INT8 pfree,
             MUMPS_INT8 * const len, MUMPS_INT8 * const iw, MUMPS_INT8 * const nv, MUMPS_INT8 * const elen,
             MUMPS_INT8 * const last, SCOTCH_Context * const esmumpscontext);
MUMPS_INT esmumpsc( const MUMPS_INT8 n, const MUMPS_INT8 iwlen, MUMPS_INT8 * const pe, const MUMPS_INT8 pfree,
             MUMPS_INT8 * const len, MUMPS_INT8 * const iw, MUMPS_INT8 * const nv, MUMPS_INT8 * const elen,
             MUMPS_INT8 * const last, SCOTCH_Context * const esmumpscontext);
#endif
#define MUMPS_SCOTCH_ORD_64 \
  F_SYMBOL(scotch_ord_64,SCOTCH_ORD_64)
void MUMPS_CALL
MUMPS_SCOTCH_ORD_64( const MUMPS_INT8 * const  n,
                 const MUMPS_INT8 * const  iwlen,
                       MUMPS_INT8 * const  petab,
                 const MUMPS_INT8 * const  pfree,
                       MUMPS_INT8 * const  lentab,
              MUMPS_INT8 * const  iwtab,
              MUMPS_INT8 * const  nvtab,
              MUMPS_INT8 * const  elentab,
              MUMPS_INT8 * const  lasttab,
              MUMPS_INT * const   ncmpa, 
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
              SCOTCH_Context * const contextptr,
#endif
              MUMPS_INT  * const  weightused,
              MUMPS_INT  * const  weightrequested );
#define MUMPS_SCOTCH_64 \
  F_SYMBOL(scotch_64,SCOTCH_64)
void MUMPS_CALL
MUMPS_SCOTCH_64( const MUMPS_INT8 * const  n,
                 const MUMPS_INT8 * const  iwlen,
                       MUMPS_INT8 * const  petab,
                 const MUMPS_INT8 * const  pfree,
                       MUMPS_INT8 * const  lentab,
                       MUMPS_INT8 * const  iwtab,
                       MUMPS_INT8 * const  nvtab,
                       MUMPS_INT8 * const  elentab,
                       MUMPS_INT8 * const  lasttab,
                       MUMPS_INT  * const  ncmpa,
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
              SCOTCH_Context * const contextptr,
#endif
                       MUMPS_INT  * const  weightused,
                       MUMPS_INT  * const  weightrequested );
#endif
#endif
