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

#include "om_format.h"
#include <ctype.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>


#define min(x,y) ((x) < (y) ? (x) : (y))

/* code by casella (https://github.com/casella), see:
 * https://github.com/OpenModelica/OpenModelica/issues/7465
 * Takes the string output of Ryu's d2s() function as input and writes
 * the minimal decimal or exponential version of it to a buffer
 * d2s_str: the string output from d2s()
 * buf: a character array (32 chars should be more than enough in all cases)
 * real_output: if true, decimal output of round numbers
 *              will have trailing .0, as in 1.0 or 1240.0
 *              if false it will be displayed as an integer 1 or 1240; additionally
 *              in case d2s returns more than 12 digits after the decimal point in
 *              the mantissa, it tries to round them to 12 digits and if there are
 *              more than 4 trailing zeros, it returns the rounded value
 */
void ryu_to_hr(const char *d2s_str, char *buf, int real_output) 
{
  int sign;  // Sign of the input
  int exp;      // Exponent in the input
  double mant;  // Mantissa in the input
  char digits[32] = {0};  // Digits of the mantissa
  int ndec = 0; // Number of digits after decimal point in the mantissa
  char str1[32] = {0}; // Stores a copy of the d2s_str input
  char str2[32] = {0}; // Stores the decimal output
  char str3[32] = {0}; // Stores the rounded decimal output, if computed
  char *ptr1 = str1;  // Movable pointer to the input copy
  char *ptr = digits; // Movable pointer to the digits string
  char *ptr2 = str2; // Movable pointer to the decimal output
  char *ptr3 = str3; // Movable pointer to the rounded mantissa
  int skip_dec = 0;  // If true, skip the decimal conversion
  int nz = 0;        // Number of trailing zeros in the rounded mantissa
  int i;  // for loop index

  // copy d2s_str input in local str1 string
  strcpy(str1, d2s_str);

  // Parse sign, digits, ndec, exp from input
  if(!strpbrk(str1, "eE")) {
    // The input string is not in (-)xx.xEyy format
    skip_dec = 1;
  } else {
    // Compute sign
    if(*ptr1 == '-') {
      sign = -1;
      ptr1++;
    } else {
      sign = +1;
    }
    // Compute digits string
    while(*ptr1 != 'E' && *ptr1 != 'e')
      *ptr++ = *ptr1++;  // Copy the digit part of the input to the digits string
    // Compute ndec
    ndec = strchr(digits, '.') ? strlen(digits) - 2 : 0;
    // Compute exp
    sscanf(++ptr1, "%d", &exp);
    // If real output with trailing .0 for round numbers is not required and
    // Ryu's mantissa has more than 12 digits after the decimal point and
    // its rounded value has at least four trailing zeros, replace Ryu's output
    // with the rounded mantissa without the trailing zeros
    if(ndec > 12 && !real_output) {
      // Compute the mantissa
      sscanf(digits, "%lf", &mant);
      // Print the rounded value of the mantissa with 12 digits after the 
      // decimal point
      sprintf(str3, "%.12f", mant);
      // Handle the case of 9.99999999999999 rounded to 10.000000000000
      if (!strcmp(str3,"10.000000000000")) {
          sprintf(str3,"1.000000000000");
          exp++;
      }
      // Remove trailing zeros from the rounded mantissa
      ptr3 += strlen(str3) - 1;
      while(*ptr3 == '0') {
        nz++;
        *ptr3-- = 0;
      }
      // Remove trailing decimal point if necessary
      if(*ptr3 == '.')
        *ptr3-- = 0;
      // Update digits string with rounded mantissa
      if(nz > 3){
        // Update digits string with rounded mantissa
        strcpy(digits, str3);
        // Update ndec
        ndec = strchr(digits, '.') ? strlen(digits) - 2 : 0;
        // Update str1 to rounded output
        ptr1 = str1;     // Rewind pointer to output copy (to update)
        ptr3 = str3;     // Rewind pointer to rounded mantissa
        if(sign == -1)
          *ptr1++ = '-'; // Write sign into str1
        while(*ptr3)
          *ptr1++ = *ptr3++;       // Copy rounded mantissa into str1
        *ptr1++ = 'e';             // Copy 'e' into str1
        sprintf(ptr1, "%d", exp);  // Print exponent into str1
      }
    }

  }
  if(exp > 5 || exp < -3) {
    // Decimal format too long, skip it (and avoid buffer overflow in str2
    skip_dec = 1;
  } else if (!skip_dec) {
    // Fill str2 with the digital form
    // Print the sign to str2
    if(sign < 0)
      *ptr2++ = '-';
    // Print the rest of the number to str2
    ptr = digits; // Rewind pointer to digits
    if(exp == 0) {
      strcpy(ptr2, digits);  // Trivial, the output is just digits
    } else if(exp > 0) { // Here you need to move the decimal point to the right
      *ptr2++ = *ptr++; // Copy first digit to str2
      ptr++; // Skip the decimal point
      for(i = 0; i < min(ndec, exp); i++)
        *ptr2++ = *ptr++; // Copy the digits up to the new decimal point
      if(exp > ndec) { // Trailing zeros are needed
        for(i = 0; i < exp - ndec; i++)
          *ptr2++ = '0'; // Add trailing zeros before decimal point
      } else if(exp < ndec) { // Remaining digits after decimal point
        *ptr2++ = '.'; // Add decimal point
        strcpy(ptr2, ptr); // Add remaining digits after decimal point
      }
    } else if(exp < 0) { // Number starts with "0.", then some zeros
      // Add "0." to str2
      *ptr2++ = '0';
      *ptr2++ = '.';
      for(i = 0; i < (-exp - 1); i++)
        *ptr2++ = '0'; // Add leading zeros after decimal point in str2
      *ptr2++ = *ptr++; // Add first digit to str2
      if(ndec > 0)
        ptr++; // Skip decimal point in digits
      strcpy(ptr2, ptr); // Add remaining digits to str2
    }
    if (exp >= ndec && real_output)  // Round decimal output, real output required
      strcat(str2, ".0");
  }
  if(exp < -3 || exp > 5 || (exp > 0 && exp - ndec > 3) || skip_dec) {
    // use the exponential form (change 'e' to 'E')
    int i = 0;
    ptr1 = str1;  // Rewind pointer to input copy
    ptr2 = buf;   // Rewind pointer to output buf
    while (*ptr1) {
      *ptr2 = *ptr1 == 'E' ? 'e' : *ptr1;
      ptr1++;
      ptr2++;
    }
    *ptr2 = '\0';  // Terminate output string
  }
  else
    strcpy(buf, str2); // use the decimal form
}

/*
 * call this one from OMEdit or other clients to print doubles (without the trailing zero)!
 * the caller needs to free the result.
 */
char* ryu_hr_tdzp(double d)
{
  char d2s_str[32] = {0};
  char buf[32] = {0};
  d2s_buffered(d, d2s_str);
  // we don't need to have the trailing zero
  ryu_to_hr(d2s_str, buf, 0);
  return strdup(buf);
}

/*
 * call this one from OMEdit or other clients to print doubles (without the trailing zero)!
 * the caller needs to provide the buffer (at least 32 chars)
 */
void ryu_hr_tdzp_buf(double d, char* buf)
{
  char d2s_str[32] = {0};
  d2s_buffered(d, d2s_str);
  // we don't need to have the trailing zero
  ryu_to_hr(d2s_str, buf, 0);
}


#if defined(TEST_RYU_TO_HR)
void test(const char *d2s, const char *expected)
{
  char buf[32];
  printf("%s -> ", d2s);
  ryu_to_hr(d2s, buf, 0);
  printf("%s", buf);
  printf("/%s\n", expected);
}

void test_real (const char *d2s, const char *expected)
{
  char buf[32];
  printf("%s -> ", d2s);
  ryu_to_hr(d2s, buf, 1);
  printf("%s", buf);
  printf("/%s\n", expected);
}

int main() 
{
  // Try all possible combinations of exp and ndec to ensure complete coverage
  printf("d2s -> output / expected\n");
  test("8e5", "8e5");
  test("8e4", "8e4");
  test("8e3", "8000");
  test("8e2", "800");
  test("8e1", "80");
  test("8e0", "8");
  test("8e-1", "0.8");
  test("8e-2", "0.08");
  test("8e-3", "0.008");
  test("8e-4", "8e-4");
  test("8e-5", "8e-5");
  test("8.13e6", "8.13e6");
  test("8.13e5", "813000");
  test("8.13e4", "81300");
  test("8.13e3", "8130");
  test("8.13e2", "813");
  test("8.13e1", "81.3");
  test("8.13e0", "8.13");
  test("8.13e-1", "0.813");
  test("8.13e-2", "0.0813");
  test("8.13e-3", "0.00813");
  test("8.13e-4", "8.13e-4");
  test("8.13e-5", "8.13e-5");
  test("8.1234567e6", "8.1234567e6");
  test("8.1234567e5", "812345.67");
  test("8.1234567e4", "81234.567");
  test("8.1234567e3", "8123.4567");
  test("8.1234567e2", "812.34567");
  test("8.1234567e1", "81.234567");
  test("8.1234567e0", "8.1234567");
  test("8.1234567e-1", "0.81234567");
  test("8.1234567e-2", "0.081234567");
  test("8.1234567e-3", "0.0081234567");
  test("8.1234567e-4", "8.1234567e-4");
  test("8.1234567e-5", "8.1234567e-5");
  test("-1.2e1", "-12");
  test("-4.56e8", "-4.56e8");
  test("1e-60", "1e-60");
  test("1e80", "1e80");
  test("NaN", "NaN");
  test("Inf", "Inf");
  test("-Inf", "-Inf");
  test("9.499999999999999e2", "950");
  test("1.9999999999999998e8", "2e8");
  test("1.9999999999999998e-6", "2e-6");
  test("-9.499999999999999e2", "-950");
  test("-1.9999999999999998e8", "-2e8");
  test("-1.9999999999999998e-6", "-2e-6");
  test("1.000000000000002e0", "1");
  test("1.0000000000000022e0", "1");
  test("-9.499999999999999e2", "-950");
  test("-1.000000000000002e0", "-1");
  test("-1.0000000000000022e0", "-1");
  test("9.99999999999999e-13", "1e-12");
  test("9.99999999999999e-5", "1e-4");
  test("9.99999999999999e-2", "0.1");
  test("9.99999999999999e-1", "1");
  test("9.99999999999999e0", "10");
  test("9.99999999999999e1", "100");
  test("9.99999999999999e5", "1e6");
  test("9.99999999999999e11", "1e12");
  test("-9.99999999999999e-13", "-1e-12");
  test("-9.99999999999999e-5", "-1e-4");
  test("-9.99999999999999e-2", "-0.1");
  test("-9.99999999999999e-1", "-1");
  test("-9.99999999999999e0", "-10");
  test("-9.99999999999999e1", "-100");
  test("-9.99999999999999e5", "-1e6");
  test("-9.99999999999999e11", "-1e12");
  test("1.234567890123456e6", "1.234567890123456e6");
  test("1.23456789012345e6", "1.23456789012345e6");
  test("1.2345678901234e6", "1.2345678901234e6");
  test("1.234567890123e6", "1.234567890123e6");
  test("1.23456789012e6", "1.23456789012e6");
  test("1.2345678901e6", "1.2345678901e6");
  test("1.23456789e6", "1.23456789e6");
  test("1.2345678e6", "1.2345678e6");
  test("1.234567e6", "1.234567e6");
  test("1.23456e6", "1.23456e6");
  test("1.2345e6", "1.2345e6");
  test("1.234e6", "1.234e6");
  test("1.23e6", "1.23e6");
  test("1.2e6", "1.2e6");
  test("1e6", "1e6");
  test("-1.234567890123456e6", "-1.234567890123456e6");
  test("-1.23456789012345e6", "-1.23456789012345e6");
  test("-1.2345678901234e6", "-1.2345678901234e6");
  test("-1.234567890123e6", "-1.234567890123e6");
  test("-1.23456789012e6", "-1.23456789012e6");
  test("-1.2345678901e6", "-1.2345678901e6");
  test("-1.23456789e6", "-1.23456789e6");
  test("-1.2345678e6", "-1.2345678e6");
  test("-1.234567e6", "-1.234567e6");
  test("-1.23456e6", "-1.23456e6");
  test("-1.2345e6", "-1.2345e6");
  test("-1.234e6", "-1.234e6");
  test("-1.23e6", "-1.23e6");
  test("-1.2e6", "-1.2e6");
  test("-1e6", "-1e6");

  test_real("8e5", "8e5");
  test_real("8e4", "8e4");
  test_real("8e3", "8000.0");
  test_real("8e2", "800.0");
  test_real("8e1", "80.0");
  test_real("8e0", "8.0");
  test_real("8e-1", "0.8");
  test_real("8e-2", "0.08");
  test_real("8e-3", "0.008");
  test_real("8e-4", "8e-4");
  test_real("8e-5", "8e-5");
  test_real("8.13e6", "8.13e6");
  test_real("8.13e5", "813000.0");
  test_real("8.13e4", "81300.0");
  test_real("8.13e3", "8130.0");
  test_real("8.13e2", "813.0");
  test_real("8.13e1", "81.3");
  test_real("8.13e0", "8.13");
  test_real("8.13e-1", "0.813");
  test_real("8.13e-2", "0.0813");
  test_real("8.13e-3", "0.00813");
  test_real("8.13e-4", "8.13e-4");
  test_real("8.13e-5", "8.13e-5");
  test_real("8.1234567e6", "8.1234567e6");
  test_real("8.1234567e5", "812345.67");
  test_real("8.1234567e4", "81234.567");
  test_real("8.1234567e3", "8123.4567");
  test_real("8.1234567e2", "812.34567");
  test_real("8.1234567e1", "81.234567");
  test_real("8.1234567e0", "8.1234567");
  test_real("8.1234567e-1", "0.81234567");
  test_real("8.1234567e-2", "0.081234567");
  test_real("8.1234567e-3", "0.0081234567");
  test_real("8.1234567e-4", "8.1234567e-4");
  test_real("8.1234567e-5", "8.1234567e-5");
  test_real("-1.2e1", "-12.0");
  test_real("-4.56e8", "-4.56e8");
  test_real("1e-60", "1e-60");
  test_real("1e80", "1e80");
  test_real("NaN", "NaN");
  test_real("Inf", "Inf");
  test_real("-Inf", "-Inf");
  test_real("9.499999999999999e2", "949.9999999999999");
  test_real("1.000000000000002e0", "1.000000000000002");
  test_real("1.0000000000000022e0", "1.0000000000000022");
  test_real("-9.499999999999999e2", "-949.9999999999999");
  test_real("-1.000000000000002e0", "-1.000000000000002");
  test_real("-1.0000000000000022e0", "-1.0000000000000022");
}
#endif // #if defined(TEST_RYU_TO_HR)
