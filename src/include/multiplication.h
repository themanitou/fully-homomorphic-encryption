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


#ifndef MULTIPLICATION_H_
#define MULTIPLICATION_H_

#include <vector>
#include <set>

using namespace std;


////////////////////////////////////////////////////
// Dang -- 2012-03-01
//
//    Data types
//
typedef set<unsigned int> set_t;
typedef vector<set_t> matrix_t;
typedef vector<matrix_t> table_t;
typedef pair<set_t::iterator, bool> set_return_t;


////////////////////////////////////////////////////
// Dang -- 2012-03-06
//
//    Global variable declarations
//
extern table_t       MultSymbTab;
extern matrix_t      MultSymbMat;
extern unsigned int  MULT_BIT_LENGTH;


////////////////////////////////////////////////////
// Dang -- 2012-03-01
//
//    Function prototypes
//
long              BuildMultSymbTable (table_t&, const unsigned int);
long              BuildMultSymbMatrix (matrix_t&, const table_t&);
long              PrintMultSymbTable (const table_t&);
long              PrintMultSymbMatrix (const matrix_t&);
long              MultSymbSets (set_t&, const set_t&, const set_t&);
long              MergeSymbSets (set_t&, const set_t&);
long              SymbCombiMN (set_t&, const matrix_t&, const unsigned int);


#endif /* MULTIPLICATION_H_ */
