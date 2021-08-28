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


#ifndef FULLHOM_H_
#define FULLHOM_H_

#include "../communication/fhe_mpi.h"

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

#include <vector>

NTL_CLIENT


////////////////////////////////////////////////////
// Dang -- 2012-02-24
//
//    Macros
//
#define max(a,b) ((a>b)?a:b)
#define min(a,b) ((a>b)?b:a)


////////////////////////////////////////////////////
// Dang -- 2012-02-23
//
//    Data types
//
typedef struct _noise_vector_t noise_vector_t, *noise_vector_pt;

struct _noise_vector_t {
    long count;
    long *elements;
};


////////////////////////////////////////////////////
// Dang -- 2012-01-31
//
//    Global variable declarations
//

// these are the domain parameters
extern ZZ          D, R, S;
extern mat_ZZ      ETA;
extern vec_ZZ      X_ARRAY;


// these variables are made global ONLY to validate different functions
extern long        LOG2_RATIO;
extern vector<long>      idxArray;
extern long        A_MAX;
extern long        SS_SIZE;
extern long        IDX_MAX;

#ifdef FHE_TEST
extern vector<long> aArray, bArray;
#endif

////////////////////////////////////////////////////
// Dang -- 2012-01-31
//
//    Function prototypes
//
void               Initialization (int, char**);
void               FHE_Initialization ();
void               CloseDown ();
void               CryptoRandomize ();
ZZX                RandPolyGen ();
ZZ                 NormalizeN (const ZZ&, const ZZ&);
long               NaifPrimaryKeyGen (ZZX&);
long               NaifKeyGen ();
long               PrimaryKeyGen (ZZX&);
long               KeyGen ();
long               CompRPowers ();
long               CompRatioPowers ();
long               PrintRPowers ();
long               PolyAtMinusX (ZZX&);
long               GModX2 (ZZ&, ZZ&, ZZX&);
noise_vector_pt    NoiseGen ();
long               FreeNoise (noise_vector_pt);
long               PrintNoise (noise_vector_pt);

long               EncryptBit (ZZ&, const long&);
long               DecryptBit (long&, const ZZ&);
long               EncryptByte (vec_ZZ&, const long&);
long               DecryptByte (long&, const vec_ZZ&);

long               DecryptWord (long&, const vec_ZZ&);

long               AndBit (ZZ&, const ZZ&, const ZZ&);
long               XorBit (ZZ&, const ZZ&, const ZZ&);
long               FlipBit (ZZ&, const ZZ&);
long               ChooseBit (ZZ&, const ZZ&, const ZZ&, const ZZ&);

long               NegByte (vec_ZZ&, const vec_ZZ&);
long               AddByte (vec_ZZ&, ZZ&, const vec_ZZ&, const vec_ZZ&);
long               SubByte (vec_ZZ&, ZZ&, const vec_ZZ&, const vec_ZZ&);
long               MaxMinSignedByte (vec_ZZ&, vec_ZZ&, const vec_ZZ&, const vec_ZZ&);
long               MulByte (vec_ZZ&, const vec_ZZ&, const vec_ZZ&);

long               AddWord (vec_ZZ&, ZZ&, const vec_ZZ&, const vec_ZZ&);

long               MulHalfByte (vec_ZZ&, const vec_ZZ&, const vec_ZZ&);
long               BuildMulBit (vec_ZZ&, const vec_ZZ&);

long               GetIndex (long&, const long&, const long&);
long               CompEtaAndXArray ();
long               RecryptBit (ZZ&, const vec_ZZ&, const vec_ZZ&, const ZZ&);
long               PostProcess (vec_ZZ&, vec_ZZ&, const ZZ&);
long               WorkerPostProcess ();
long               WorkerRecrypt ();

#ifdef FHE_DEBUG_LOG
long               PrintAllVariables ();
#endif

#ifdef FHE_TEST
void               FHE_Partial_Initialization ();
long               ReCompXArray ();
long               ReCompEta ();
#endif

#endif /* FULLHOM_H_ */
