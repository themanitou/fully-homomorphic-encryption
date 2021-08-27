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


#include "multiplication.h"
#include "fullhom.h"


////////////////////////////////////////////////////
// Dang -- 2012-03-06
//
//    Global variable creations
//
table_t        MultSymbTab;
matrix_t       MultSymbMat;
unsigned int   MULT_BIT_LENGTH = 4;


//////////////////////////////////////////////////////
// Dang -- 2012-03-05
// Function description:
//       tbc...
//
long MergeSymbSets (set_t& result, const set_t& aSet) {

  set_t::iterator i;
  set_return_t ret;

  for (i = aSet.begin (); i != aSet.end (); i++) {
    ret = result.insert (*i);

    if (ret.second == false) { // if this element already exists then remove it
      result.erase (ret.first);
    }
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-03-05
// Function description:
//       tbc...
//
long SymbCombiMN (set_t& result, const matrix_t& aColumn,
		  const unsigned int m) {

  unsigned int n;
  set_t aSet;

  vector<unsigned int> aVector;
  unsigned int i, j, k;
  unsigned int count;

  n = aColumn.size ();
  if (n < m) {
    cout << "[SymbCombiMN] Error: n < m." << endl;
    return -1;
  }

  // initialize the combination vector
  // aVector[0] = 1, aVector[1] = 2, ..., aVector[m-1] = m
  for (i = 0; i < m; ) {
    aVector.push_back (++i);
  }

  count = 0;

  // iterate until we have got all combinations
  while (1) {

    // print the current combination
    count++;
    aSet = aColumn[aVector[0] - 1];
    for (i = 1; i < m; i++) {
      MultSymbSets (aSet, aSet, aColumn[aVector[i] - 1]);
    }
    // merge the combination set into the result set
    MergeSymbSets (result, aSet);
    // clear the combination set
    aSet.clear ();

    // check if we have enumerated all combinations
    if (aVector[0] == (n - m + 1)) {
      break;
    }

    // start looking from the tail to the head
    for (i = 0; i < m; i++) {
      // k is the tail
      k = m - 1 - i;
      // if we can still increase by 1
      if (aVector[k] < (n - i)) {
	// then increase by 1
	aVector[k]++;
	// and reset the rest of the tail to the next element just like in
	// the above initialization, except that we shift by some amount
	for (j = k + 1; j < m; j++) {
	  aVector[j] = aVector[j - 1] + 1;
	}
	break; // break out of for (i = 0; ...)
      }
    } // end for (i = 0, ...)

  } // end while (1)

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-03-05
// Function description:
//       tbc...
//
long MultSymbSets (set_t& result, const set_t& a, const set_t& b) {

  unsigned int c;
  set_t aSet;
  set_t::iterator i, j;
  set_return_t ret;

  for (i = a.begin (); i != a.end (); i++) {
    for (j = b.begin (); j != b.end (); j++) {

      c = ((*i) | (*j));
      ret = aSet.insert (c);
      if (ret.second == false) { // this symbol already exists, then remove it
	aSet.erase (ret.first);
      }

    }
  }

  result = aSet;

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-03-01
// Function description:
//       tbc...
//
long PrintMultSymbTable (const table_t& symbTab) {

  matrix_t col;
  set_t aSet;

  unsigned int i, j;
  set_t::iterator k;

  cout << hex;
  for (i = 0; i < symbTab.size (); i++) {
    cout << "[PrintMultSymbTable] Column #" << i << " : ";
    col = symbTab[i];

    for (j = 0; j < col.size (); j++) {
      cout << " + [";
      aSet = col[j];

      for (k = aSet.begin (); k != aSet.end (); k++) {
	cout << "(" << (*k) << "), ";
      }
      cout << "]";
    }
    cout << endl;
  }
  cout << dec;

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-03-01
// Function description:
//       tbc...
//
long BuildMultSymbTable (table_t& symbTab, const unsigned int bitLen) {

  matrix_t col;
  set_t aSet;

  unsigned int maxLen;
  unsigned int i, j, m, n;
  unsigned int lowBit, highBit;

  symbTab.clear ();

  maxLen = (bitLen << 1) - 1;
  lowBit = 1;
  highBit = (1 << bitLen);

  // build the first-half base symbols
  for (i = 0; i < bitLen; i++) {
    for (j = 0; j <= i; j++) {
      // b_{j}a_{i-j}
      aSet.insert ((highBit << j) | (lowBit << (i - j)));
      col.push_back (aSet);
      aSet.clear ();
    }
    symbTab.push_back (col);
    col.clear ();
  }

  // build the second-half base symbols
  for (i = bitLen; i < maxLen; i++) {
    for (j = (i - bitLen) + 1; j < bitLen; j++) {
      // b_{j}a_{i-j}
      aSet.insert ((highBit << j) | (lowBit << (i - j)));
      col.push_back (aSet);
      aSet.clear ();
    }
    symbTab.push_back (col);
    col.clear ();
  }

  symbTab.resize (symbTab.size () + 1);

  // build the carry bit symbols
  for (i = 0; i < maxLen; i++) {
    col = symbTab[i];
    n = col.size ();
    for (j = i + 1; j <= maxLen; j++) {
      m = (1 << (j - i));
      if (n < m) {
	break;
      }

      SymbCombiMN (aSet, col, m);
      symbTab[j].push_back (aSet);
      aSet.clear ();
    }
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-03-06
// Function description:
//       tbc...
//
long BuildMultSymbMatrix (matrix_t& symbMat, const table_t& symbTab) {

  unsigned int i, j;
  set_t aSet;
  matrix_t col;

  symbMat.clear ();

  for (i = 0; i < symbTab.size (); i++) {
    col = symbTab[i];

    for (j = 0; j < col.size (); j++) {
      MergeSymbSets (aSet, col[j]);
    }

    symbMat.push_back (aSet);
    aSet.clear ();
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-03-06
// Function description:
//       tbc...
//
long PrintMultSymbMatrix (const matrix_t& symbMat) {

  set_t aSet;

  unsigned int i;
  set_t::iterator k;

  cout << hex;
  for (i = 0; i < symbMat.size (); i++) {

    aSet = symbMat[i];
    cout << "[PrintMultSymbMatrix] Column #" << i << " : [";

    for (k = aSet.begin (); k != aSet.end (); k++) {
      cout << "(" << (*k) << "), ";
    }

    cout << "]" << endl;
  }
  cout << dec;

  return 0;
}
