/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF THE BSD NEW LICENSE OR THE
 * GPL VERSION 3 LICENSE OR THE OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the OSMC (Open Source Modelica Consortium)
 * Public License (OSMC-PL) are obtained from OSMC, either from the above
 * address, from the URLs: http://www.openmodelica.org or
 * http://www.ida.liu.se/projects/OpenModelica, and in the OpenModelica
 * distribution. GNU version 3 is obtained from:
 * http://www.gnu.org/copyleft/gpl.html. The New BSD License is obtained from:
 * http://www.opensource.org/licenses/BSD-3-Clause.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, EXCEPT AS
 * EXPRESSLY SET FORTH IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE
 * CONDITIONS OF OSMC-PL.
 *
 */


#ifndef OM_FORMAT_H
#define OM_FORMAT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ryu.h"

/* code by casella (https://github.com/casella), see:
 * https://github.com/OpenModelica/OpenModelica/issues/7465
 * Takes the string output of Ryu's d2s() function as input and writes
 * the minimal decimal or exponential version of it to a buffer
 * str1: the string output from d2s()
 * buf: a character array (32 chars should be more than enough in all cases)
 * real_output: if true, decimal output of round numbers
 *              will have trailing .0, as in 1.0 or 1240.0
 *              if false it will be displayed as an integer 1 or 1240
 */
extern void ryu_to_hr(const char *str1, char *buf, int real_output);

/*
 * call this one from OMEdit or other clients to print doubles (without the trailing zero)!
 * the caller needs to free the result.
 */
extern char* ryu_hr_tdzp(double r);

/*
 * call this one from OMEdit or other clients to print doubles (without the trailing zero)!
 * the caller needs to provide the buffer (at least 32 chars)
 */
extern void ryu_hr_tdzp_buf(double r, char* buf);

#ifdef __cplusplus
}
#endif

#endif // OM_FORMAT_H
