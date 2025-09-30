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
/* $Id */
#ifndef MUMPS_FLYTES_H
#define MUMPS_FLYTES_H
#include <stdint.h>
#include "mumps_common.h"
#if !defined(USE_AVX512_VBMI)
#undef __AVX512F__
#undef __AVX512VBMI__
#endif
void MUMPS_CALL mumps_flyte_return();
#endif
