/******************************************************************
 *
 *   Fully-Homomorphic Cryptography library,
 *   based on Gentry-Halevi ideal lattice scheme.
 *
 *   Author: Quan Nguyen (https://github.com/themanitou)
 *
 *   This library is open-source software distributed under the
 *   terms of the GNU Lesser General Public License (LGPL) version
 *   2.1 or later.  See the file doc/copying.txt for complete
 *   details on the licensing of this library.
 *
 *******************************************************************/


#include "../include/tools.h"


//////////////////////////////////////////////////////
// Dang -- 2012-05-28
// Function description:
//       tbc...
//
string convBase (const unsigned long& x, const long& base) {

  string digits = "0123456789abcdef";
  string result;
  unsigned long v;

  v = x;

  if((base < 2) || (base > 16))
    result = "[convBase] Error: base out of range.";
  else
    do {
      result = digits[v % base] + result;
      v /= base;
    }
    while(v);

  return result;
}
