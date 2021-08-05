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


#include <sys/time.h>
#include "../include/fhe_mpi.h"
#include "../include/tools.h"
#include "../include/multiplication.h"
#include "../include/fullhom.h"


////////////////////////////////////////////////////
// Dang -- 2012-01-31
//
//    Global variable creations
//

// These variables are declared EXTERN since they are the domain parameters
ZZ          D, R, S;  // PK = (D,R), SK = S
mat_ZZ      ETA;                // matrix of eta values, only two of each row
                                // are encryption of 1, the others are enc of 0
vec_ZZ      X_ARRAY;            // representative elements of the big sets

// These variables are NOT declared EXTERN
long        LOG2_N;
long        N;
long        MAX_BIT_LENTH; // t-bit length coefficients
ZZ          MAX_RANDOM;

ZZX         F;  // F = X^N + 1
long        NOISE_MAGNITUDE_MIN, NOISE_MAGNITUDE_MAX;

long        SS_SIZE;            // sparse subset size -- s
long        BS_SIZE;            // big set size -- S
long        IDX_MAX;            // maximum index -- c
long        LOG2_RATIO;         // log(ratio between elmts in big set)
long        PREC_PARAM;         // preceision parameter -- p

vec_ZZ      R_POWERS;
long        A_MAX;
ZZ          RATIO;
vec_ZZ      RATIO_POWERS;
RR          D_PREC;

// this variable is global ONLY on the purpose OF validation
// it should be local within the CompEtaAndXArray function
vector<long> idxArray;  // this array contain the indexes i of the chunk X_k^i
                        // in the big set B_k; k = 1 .. SS_SIZE

#ifdef FHE_TEST
vector<long> aArray, bArray;
#endif


#ifdef FHE_DEBUG_LOG
//////////////////////////////////////////////////////
// Dang -- 2012-06-13
// Function description:
//       tbc...
//
long PrintAllVariables () {
  // print D
  cout << "D=" << D << endl;

  // print R
  cout << "R=" << R << endl;

  // print S
  cout << "S=" << S << endl;

  // print X_ARRAY
  cout << "X_ARRAY=" << X_ARRAY << endl;

  // print ETA
  cout << "ETA=" << ETA << endl;

  return 0;
}
#endif

//////////////////////////////////////////////////////
// Dang -- 2012-06-07
// Function description:
//       tbc...
//
long BcastParams () {

  long i, j;

  // broadcast D
  MPI_ZZ_Bcast (D);
#ifdef FHE_MPI_LOG
  cout << "[BcastParams-" << FHE_MPI_RANK << "] Broadcast D=" << trunc_long (D, 32) << endl;
#endif

  // compute D_PREC
  if (FHE_MPI_RANK)
    conv (D_PREC, D);

  // broadcast R
  MPI_ZZ_Bcast (R);
#ifdef FHE_MPI_LOG
  cout << "[BcastParams-" << FHE_MPI_RANK << "] Broadcast R=" << trunc_long (R, 32) << endl;
#endif

  // slaves compute R_POWERS and RATIO_POWERS
  if (FHE_MPI_RANK) {
    CompRPowers ();
    CompRatioPowers ();
  }

  // begin test X_ARRAY
  // broadcast S
  MPI_ZZ_Bcast (S);
#ifdef FHE_MPI_LOG
  cout << "[BcastParams-" << FHE_MPI_RANK << "] Broadcast S=" << trunc_long (S, 32) << endl;
#endif

  ZZ z;
  if (FHE_MPI_RANK)
    idxArray.resize (SS_SIZE);
  for (i = 1; i <= SS_SIZE; i++) {
    z = idxArray[i - 1];
    MPI_ZZ_Bcast (z);
    idxArray[i - 1] = trunc_long (z, 32);
  }
#ifdef FHE_MPI_LOG
  cout << "[BcastParams-" << FHE_MPI_RANK << "] Broadcast idxArray." << endl;
#endif
  // end test X_ARRAY

  // broadcast X_ARRAY
  for (i = 1; i <= SS_SIZE; i++)
    MPI_ZZ_Bcast (X_ARRAY(i));
#ifdef FHE_MPI_LOG
  cout << "[BcastParams-" << FHE_MPI_RANK << "] Broadcast X_ARRAY." << endl;
#endif

  // broadcast ETA
  for (i = 1; i <= SS_SIZE; i++)
    for (j = 1; j <= IDX_MAX; j++)
      MPI_ZZ_Bcast (ETA(i, j));
#ifdef FHE_MPI_LOG
  cout << "[BcastParams-" << FHE_MPI_RANK << "] Broadcast ETA." << endl;
#endif

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
void Initialization (int argc, char* argv[]) {

  // MPI initialization
//  MPI_Init (&argc, &argv);
//  MPI_Comm_size(MPI_COMM_WORLD, &FHE_MPI_SIZE);
//  MPI_Comm_rank(MPI_COMM_WORLD, &FHE_MPI_RANK);
  MPI::Init (argc, argv);
  FHE_MPI_SIZE = MPI::COMM_WORLD.Get_size ();
  FHE_MPI_RANK = MPI::COMM_WORLD.Get_rank ();

  // random seed
  CryptoRandomize ();

  // set initial values
  LOG2_N = 9;
  N = (1 << LOG2_N);
  MAX_BIT_LENTH = 380; // t-bit length coefficients
  MAX_RANDOM = power2_ZZ (MAX_BIT_LENTH);
  F = ZZX (N, 1) + ZZX (0, 1);  // F = X^N + 1

  NOISE_MAGNITUDE_MIN = 15;  // based on the fact that N = 128,
  NOISE_MAGNITUDE_MAX = 20;  // NOISE_MAGNITUDE should
                             // change if N is different

  SS_SIZE    = 15;    // sparse subset size -- s
  BS_SIZE    = 512;   // big set size -- S
  IDX_MAX    = ceil (2.0 * sqrt (BS_SIZE)); // maximum index -- c
  LOG2_RATIO = 26;    // log(ratio between elmts in big set)
  PREC_PARAM = 4;     // preceision parameter -- p

  A_MAX = ceil (((2.0*IDX_MAX - 1) - sqrt ((2.0*IDX_MAX - 1)*(2.0*IDX_MAX - 1) - 8.0*BS_SIZE)) / 2.0);
  RATIO = power2_ZZ (LOG2_RATIO);

  // allocate memory for vectors and matrices
  R_POWERS.SetLength (N);
  X_ARRAY.SetLength (SS_SIZE);
  ETA.SetDims (SS_SIZE, IDX_MAX);
  RATIO_POWERS.SetLength (BS_SIZE);

  // FHE initialization
  FHE_Initialization ();

}


#ifdef FHE_TEST
//////////////////////////////////////////////////////
// Dang -- 2012-06-13
// Function description:
//       tbc...
//
void FHE_Partial_Initialization () {

  // only the master process compute the domain parameters
  if (FHE_MPI_RANK == 0) {
    // compute ETA matrix and X array
    CompEtaAndXArray ();
  }

  // broadcast the domain parameters
  BcastParams ();

}
#endif


//////////////////////////////////////////////////////
// Dang -- 2012-06-07
// Function description:
//       tbc...
//
void FHE_Initialization () {

  // only the master process compute the domain parameters
  if (FHE_MPI_RANK == 0) {

    // print the irreducible polynomial
#ifdef FHE_LOG
    cout << "[Master.Initialization] F=" << F << endl;
#endif

    // generate keys
    KeyGen ();

    // compute R_POWERS
    CompRPowers ();

    // compute the ratio powers
    CompRatioPowers ();

#ifdef FHE_MULT_SYMB
    // compute multiplicative table of symbols for MULT_BIT_LENGTH bits
    BuildMultSymbTable (MultSymbTab, MULT_BIT_LENGTH);

#ifdef FHE_LOG
    PrintMultSymbTable (MultSymbTab);
#endif

    BuildMultSymbMatrix (MultSymbMat, MultSymbTab);

#ifdef FHE_LOG
    PrintMultSymbMatrix (MultSymbMat);
#endif
#endif

    // compute ETA matrix and X array
    CompEtaAndXArray ();

  }

  // broadcast the domain parameters
  BcastParams ();

}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
void CloseDown () {
  // MPI shut down
  MPI_Finalize ();
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long CompRPowers () {

  long i;

  R_POWERS(1) = 1;
  R_POWERS(2) = R;
  for (i = 3; i <= N; i++) {
    R_POWERS(i) = MulMod (R, R_POWERS(i-1), D);
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-08
// Function description:
//       tbc...
//
long CompRatioPowers () {

  long i;

  RATIO_POWERS(1) = RATIO;
  for (i = 2; i <= BS_SIZE; i++) {
    MulMod (RATIO_POWERS(i), RATIO, RATIO_POWERS(i-1), D);
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long PrintRPowers () {

  long i;

  cout << "[PrintRPowers] R_POWERS=[";
  for (i = 1; i <= N; i++)
    cout << R_POWERS(i) << ",";
  cout << "]." << endl;

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-01-31
// Function description:
//       tbc...
//
void CryptoRandomize () {

  ZZ x;

  struct timeval tv;
  gettimeofday (&tv, NULL);

  ZZFromBytes (x, (unsigned char *)&(tv.tv_usec), sizeof (tv.tv_usec));

  SetSeed (x);
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
noise_vector_pt NoiseGen () {

  noise_vector_pt pNoise;
  long i, j;
  ZZ x;

  pNoise = (noise_vector_pt)malloc (sizeof (noise_vector_t));

  if (pNoise == NULL) {
    // Error, cannot allocate noise vector
    cout << "[NoiseGen] Error, cannot allocate noise vector." << endl;
    return NULL;
  }

  pNoise->count = NOISE_MAGNITUDE_MIN + RandomBnd (NOISE_MAGNITUDE_MAX - NOISE_MAGNITUDE_MIN + 1);
  pNoise->elements = (long *)malloc (pNoise->count * sizeof (long));

  if (pNoise->elements == NULL) {
    // Error, cannot allocate noise elements
    cout << "[NoiseGen] Error, cannot allocate noise elements." << endl;
    return NULL;
  }

  x = 0;

  for (i = 0; i < pNoise->count; i++) {
    do {
      j = RandomBnd (N);
    } while (bit (x, j) == 1);
    SetBit (x, j);

    if (RandomBnd (2) == 1)
      j = -j;

    pNoise->elements[i] = j;
  }

  return pNoise;

}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long FreeNoise (noise_vector_pt pNoise) {

  if (pNoise) {

    if (pNoise->elements && pNoise->count)
      free (pNoise->elements);

    free (pNoise);
    pNoise = NULL;
  }

  return 0;

}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long PrintNoise (noise_vector_pt pNoise) {
  long i;

  if (pNoise) {
    cout << "[PrintNoise] count=" << pNoise->count;
    if (pNoise->count) {
      cout << "  [";
      for (i = 0; i < pNoise->count; i++)
	cout << pNoise->elements[i] << ",";
      cout << "]";
    }
    cout << "." << endl;
  }
  else {
    cout << "[PrintNoise] Noise vector is not created." << endl;
    return 1;
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-01-31
// Function description:
//       tbc...
//
ZZX RandPolyGen () {

  long i;
  ZZX result;

  for (i = 0; i < N; i++)
    SetCoeff (result, i, RandomBnd (MAX_RANDOM) - (MAX_RANDOM >> 1));

  result.normalize ();

  return result;
}


//////////////////////////////////////////////////////
// Dang -- 2012-01-31
// Function description:
//       tbc...
//
ZZ NormalizeN (const ZZ& zzA, const ZZ& zzN) {

  if ((zzA << 1) >= zzN)
    return zzA - zzN;
  else
    return zzA;
}


//////////////////////////////////////////////////////
// Dang -- 2012-01-31
// Function description:
//       tbc...
//
long NaifPrimaryKeyGen (ZZX& v) {

  ZZX w, t;
  ZZ  w0_inv, r_pow_n, w0, w1;
  long inv_status, i;
  bool secretKeyFound = false;

  // compute the scaled inverse of V, namely W such that V*W = D mod F
  XGCD (D, w, t, v, F, 1);

  // if D is even then regenerate
  if (!IsOdd (D)) {
    return 1;
  }

  // initialize D_PREC
  conv (D_PREC, D);

  w0 = ConstTerm (w);
  w1 = coeff (w, 1);

  // compute R
  if (w0 < 0)
    inv_status = InvModStatus (w0_inv, w0 + D, D);
  else
    inv_status = InvModStatus (w0_inv, w0, D);

  if (inv_status) {
    return 2;
  }

  MulMod (R, w1, w0_inv, D);
  r_pow_n = NormalizeN (PowerMod (R, N, D), D);

  if (r_pow_n != -1) {
    return 3;
  }

  // search for an odd coefficient w_i to be the secret key
  i = 0;
  while (i <= deg (w)) {
    if (IsOdd (coeff (w, i))) {
      S = coeff (w, i);
      secretKeyFound = true;
      break;
    }
    i++;
  }

  if (!secretKeyFound) {
    return 4;
  }

  return 0;

}


//////////////////////////////////////////////////////
// Dang -- 2012-02-02
// Function description:
//       tbc...
//
long NaifKeyGen () {

  ZZX v;
  long i;

  do {
    // Generate v
    v = RandPolyGen ();

    i = NaifPrimaryKeyGen (v);

    switch (i) {

    case 1:
      cout << "[NaifKeyGen] D is even.\n";
      break;

    case 2:
      cout << "[NaifKeyGen] w0 does not have inverse mod D.\n";
      break;

    case 3:
      cout << "[NaifKeyGen] R^N != -1 Mod D.\n";
      break;

    case 4:
      cout << "[NaifKeyGen] Cannot find secret key.\n";
      break;

    case 0:
      cout << "[NaifKeyGen] Keys generated successfully.\n";
      break;

    default:
      break;

    }

  } while (i);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-01-31
// Function description:
//       tbc...
//
long PolyAtMinusX (ZZX& v) {

  long i;

  for (i = 1; i <= deg (v); i = i + 2)
    SetCoeff (v, i, -coeff (v, i));

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-01-31
// Function description:
//       tbc...
//
long GModX2 (ZZ& g0, ZZ& g1, ZZX& V) {

  ZZX u, u_mx, v, v_mx, a, b, f;
  long n_j, i;

  v = ZZX (V);
  u = ZZX (0, 1);
  n_j = N;

  while (n_j >= 2) {

    f = ZZX (n_j, 1) + ZZX (0, 1);

    u_mx = u;
    v_mx = v;
    PolyAtMinusX (u_mx);
    PolyAtMinusX (v_mx);

    a = MulMod (v, v_mx, f);
    b = MulMod (u, v_mx, f) + MulMod (u_mx, v, f);

    n_j >>= 1;
    clear (u);
    clear (v);

    for (i = 0; i < n_j; i++) {
      SetCoeff (u, i, coeff (b, i << 1));
      SetCoeff (v, i, coeff (a, i << 1));
    }
  }

  // come to this stage, degree of u, v should be zero
  if (deg (v) || deg (u)) {
    cout << "[GModX2] Error, degree u=" << deg(u) << ", degree v=" << deg (v) << "\n";
    return 1;
  }

  g0 = ConstTerm (v);
  g1 = ConstTerm (u);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-01-31
// Function description:
//       tbc...
//
long PrimaryKeyGen (ZZX& v) {

  ZZ  w0, w1, t, w0_inv, r_pow_n;
  long i, inv_status;
  bool secretKeyFound = false;

  i = GModX2 (D, w0, v);
  if (i) {
    return -1;
  }

  if (!IsOdd (D)) {
    return 1;
  }

  // initialize D_PREC
  conv (D_PREC, D);  

  w0 >>= LOG2_N;   // w0 = g1 / N;

  v = MulByXMod (v, F);  // v = X*v mod F

  i = GModX2 (t, w1, v);
  if (i) {
    return -2;
  }

  w1 >>= LOG2_N;   // w1 = g1 / N;

  // compute R
  if (w0 < 0)
    inv_status = InvModStatus (w0_inv, w0 + D, D);
  else
    inv_status = InvModStatus (w0_inv, w0, D);

  if (inv_status) {
    return 2;
  }

  MulMod (R, w1, w0_inv, D);
  r_pow_n = NormalizeN (PowerMod (R, N, D), D);

  if (r_pow_n != -1) {
    return 3;
  }

  // search for an odd coefficient w_i to be the secret key
  if (IsOdd (w0)) {
    if (w0 < 0)
      S = w0 + D;
    else
      S = w0;
    secretKeyFound = true;
  }
  else if (IsOdd (w1)) {
    if (w1 < 0)
      S = w1 + D;
    else
      S = w1;
    secretKeyFound = true;
  }
  else {
    i = 2;
    t = w1;
    while (i <= N) {
      //      t = NormalizeN (MulMod (t, R, D), D);
      MulMod (t, t, R, D);
      if (IsOdd (NormalizeN (t, D))) {
	S = t;
	secretKeyFound = true;
	break;
      }
      i++;
    }
  }

  if (!secretKeyFound) {
    return 4;
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-02
// Function description:
//       tbc...
//
long KeyGen () {

  ZZX v;
  long i;

  do {
    // Generate v
    v = RandPolyGen ();

    i = PrimaryKeyGen (v);

    switch (i) {

    case -1:
      cout << "[KeyGen] GmodX2 1st call's error.\n";
      break;
    
    case -2:
      cout << "[KeyGen] GmodX2 2nd call's error.\n";
      break;

    case 1:
      cout << "[KeyGen] D is even.\n";
      break;

    case 2:
      cout << "[KeyGen] w0 does not have inverse mod D.\n";
      break;

    case 3:
      cout << "[KeyGen] R^N != -1 Mod D.\n";
      break;

    case 4:
      cout << "[KeyGen] Cannot find secret key.\n";
      break;

    case 0:
      cout << "[KeyGen] Keys generated successfully.\n";
      break;

    default:
      break;

    }

  } while (i);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long EncryptBit (ZZ& result, const long& aBit) {

  noise_vector_pt pNoise;
  long i, j;

  pNoise = NoiseGen ();

  result = 0;
  for (i = 0; i < pNoise->count; i++) {
    j = pNoise->elements[i];
    if (j > 0) {
      result += R_POWERS[j];
    }
    else if (j < 0) {
      result -= R_POWERS[-j];
    }
    if (result >= D)
      result -= D;
    else if (result < 0)
      result += D;
  }
  result = ((result << 1) + (aBit & 0x1)) % D;

  //  result = NormalizeN (result, D);

  FreeNoise (pNoise);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long DecryptBit (long& result, const ZZ& aCiphertext) {

  ZZ temp;

  temp = aCiphertext;
  if (aCiphertext < 0)
    temp += D;

  if (temp < 0 || temp >= D) {
    // Error, invalid ciphertext
    cout << "[DecryptBit] Error, invalid ciphertext." << endl;
    return 1;
  }

  temp = NormalizeN (MulMod (temp, S, D), D);

  result = IsOdd (temp);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long EncryptByte (vec_ZZ& result, const long& aByte) {

  long i;

  result.SetLength (8);

  for (i = 0; i < 8; i++) {
    EncryptBit (result[i], aByte >> i);
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long DecryptByte (long& result, const vec_ZZ& ciphertext) {

  long i, temp;

  if (ciphertext.length () != 8) {
    // Error, ciphertext invalid
    cout << "[DecryptByte] Error, invalid ciphertext." << endl;
    return 1;
  }

  result = 0;
  for (i = 0; i < 8; i++) {
    DecryptBit (temp, ciphertext[i]);
    result |= (temp << i);
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long AndBit (ZZ& result, const ZZ& a, const ZZ& b) {

  result = MulMod (a, b, D);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-03-08
// Function description:
//       tbc...
//
long XorBit (ZZ& result, const ZZ& a, const ZZ& b) {

  result = AddMod (a, b, D);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-24
// Function description:
//       result = xa + (1-x)b
//
long ChooseBit (ZZ& result, const ZZ& x, const ZZ& a, const ZZ& b) {

  ZZ x_flip, t1, t2;

  FlipBit (x_flip, x);
  AndBit (t1, x, a);
  AndBit (t2, x_flip, b);
  XorBit (result, t1, t2);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long FlipBit (ZZ& result, const ZZ& a) {

  result = AddMod (a, 1, D);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long NegByte (vec_ZZ& result, const vec_ZZ& a) {

  long i;
  ZZ c;
  vec_ZZ one;

  if (a.length () != 8) {
    // Error, ciphertext invalid
    cout << "[NegByte] Error, invalid ciphertext." << endl;
    return 1;
  }

  result.SetLength (8);

  for (i = 0; i < 8; i++) {
    FlipBit (result[i], a[i]);
  }

  one.SetLength (8);
  clear (one);
  one[0] = 1;
  AddByte (result, c, result, one);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long AddByte (vec_ZZ& result, ZZ& carry, const vec_ZZ& a, const vec_ZZ& b) {

  long i;
  ZZ c, t0, t1, t2, t3;

  if (a.length () != 8 || b.length () != 8) {
    // Error, ciphertext invalid
    cout << "[AddByte] Error, invalid ciphertext." << endl;
    return 1;
  }

  c = 0;
  result.SetLength (8);

  for (i = 0; i < 8; i++) {
    XorBit (t0, a[i], b[i]);
    XorBit (t0, t0, c);

    AndBit (t1, a[i], b[i]);
    AndBit (t2, a[i], c);
    AndBit (t3, b[i], c);

    result[i] = t0;

    XorBit (c, t1, t2);
    XorBit (c, c, t3);
  }

  carry = c;

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-24
// Function description:
//       tbc...
//
long SubByte (vec_ZZ& result, ZZ& carry, const vec_ZZ& a, const vec_ZZ& b) {

  vec_ZZ c;

  if (a.length () != 8 || b.length () != 8) {
    // Error, ciphertext invalid
    cout << "[SubByte] Error, invalid ciphertext." << endl;
    return 1;
  }

  result.SetLength (8);

  NegByte (c, b);
  AddByte (result, carry, a, c);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-24
// Function description:
//       tbc...
//
long MaxMinSignedByte (vec_ZZ& bMax, vec_ZZ& bMin, const vec_ZZ& a, const vec_ZZ& b) {

  vec_ZZ c;
  ZZ carry;
  ZZ t0, t1, t2, t3, t4;
  long i;

  if (a.length () != 8 || b.length () != 8) {
    // Error, ciphertext invalid
    cout << "[MaxMinByte] Error, invalid ciphertext." << endl;
    return 1;
  }

  XorBit (t0, a[7], b[7]);

  c.SetLength (8);
  SubByte (c, carry, a, b);

  bMax.SetLength (8);
  bMin.SetLength (8);

  for (i = 0; i < 8; i++) {
    ChooseBit (t1, a[7], b[i], a[i]);
    ChooseBit (t2, c[7], b[i], a[i]);
    ChooseBit (t3, t0, t1, t2);

    ChooseBit (t1, a[7], a[i], b[i]);
    ChooseBit (t2, c[7], a[i], b[i]);
    ChooseBit (t4, t0, t1, t2);

    bMax[i] = t3;
    bMin[i] = t4;
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-24
// Function description:
//       tbc...
//
long AddWord (vec_ZZ& result, ZZ& carry, const vec_ZZ& a, const vec_ZZ& b) {

  long i;
  ZZ c, t0, t1, t2, t3;

  if (a.length () != 16 || b.length () != 16) {
    // Error, ciphertext invalid
    cout << "[AddWord] Error, invalid ciphertext." << endl;
    return 1;
  }

  c = 0;
  result.SetLength (16);

  for (i = 0; i < 16; i++) {
    XorBit (t0, a[i], b[i]);
    XorBit (t0, t0, c);

    AndBit (t1, a[i], b[i]);
    AndBit (t2, a[i], c);
    AndBit (t3, b[i], c);

    result[i] = t0;

    XorBit (c, t1, t2);
    XorBit (c, c, t3);
  }

  carry = c;

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-24
// Function description:
//       tbc...
//
long MulByte (vec_ZZ& result, const vec_ZZ& a, const vec_ZZ& b) {

  vec_ZZ c;
  ZZ carry;
  long i, j;

  if (a.length () != 8 || b.length () != 8) {
    // Error, ciphertext invalid
    cout << "[MulByte] Error, invalid ciphertext." << endl;
    return 1;
  }

  c.SetLength (16);
  result.SetLength (16);
  clear (result);

  for (i = 0; i < 4; i++) {

    clear (c);

    for (j = 0; j < 4; j++)
      AndBit (c[i+j], b[i], a[j]);

    AddWord (result, carry, result, c);
    result[i + 8] = carry;
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-24
// Function description:
//       tbc...
//
long DecryptWord (long& result, const vec_ZZ& ciphertext) {

  long i, temp;

  if (ciphertext.length () != 16) {
    // Error, ciphertext invalid
    cout << "[DecryptWord] Error, invalid ciphertext." << endl;
    return 1;
  }

  result = 0;
  for (i = 0; i < 16; i++) {
    DecryptBit (temp, ciphertext[i]);
    result |= (temp << i);
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-03-06
// Function description:
//       tbc...
//
long BuildMulBit (vec_ZZ& result, const vec_ZZ& a) {

  if (a.length () != 8) {
    // Error, ciphertext invalid
    cout << "[BuildMulBit] Error, invalid ciphertext." << endl;
    return 1;
  }

  result.SetLength (16);
  clear (result);

  result[1] = a[0];
  result[2] = a[1];
  result[4] = a[2];
  result[8] = a[3];

  AndBit (result[3], result[1], result[2]);
  AndBit (result[5], result[1], result[4]);
  AndBit (result[6], result[2], result[4]);
  AndBit (result[7], result[2], result[5]);

  AndBit (result[9], result[1], result[8]);
  AndBit (result[10], result[2], result[8]);
  AndBit (result[11], result[3], result[8]);
  AndBit (result[12], result[4], result[8]);

  AndBit (result[13], result[5], result[8]);
  AndBit (result[14], result[6], result[8]);
  AndBit (result[15], result[7], result[8]);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-03-06
// Function description:
//       tbc...
//
long MulHalfByte (vec_ZZ& result, const vec_ZZ& a, const vec_ZZ& b) {

  ZZ c, d;
  vec_ZZ mulBitA, mulBitB;

  unsigned int i, highIdx, lowIdx;
  set_t aSet;
  set_t::iterator j;

  if (a.length () != 8 || b.length () != 8) {
    // Error, ciphertext invalid
    cout << "[MulHalfByte] Error, invalid ciphertext." << endl;
    return 1;
  }

  result.SetLength (8);
  clear (result);

  BuildMulBit (mulBitA, a);
  BuildMulBit (mulBitB, b);

  for (i = 0; i < MultSymbMat.size (); i++) {
    aSet = MultSymbMat[i];
    clear (c);

    for (j = aSet.begin (); j != aSet.end (); j++) {
      clear (d);

      highIdx = (*j) >> 4;
      lowIdx = (*j) & 0xF;

      AndBit (d, mulBitB[highIdx], mulBitA[lowIdx]);
      XorBit (c, c, d);
    }

    result[i] = c;
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-05-22
// Function description:
//       tbc...
//
long GetIndex (long& result, const long& a, const long& b) {

  if ((a < 1) || (a > A_MAX) ||
      (b < 1) || (b > IDX_MAX) ||
      (a == b))
    return -1;
  else
    result = (a - 1)*IDX_MAX - ((a*(a - 1)) >> 1) + (b - a);

  return 0;
}


#ifdef FHE_TEST
//////////////////////////////////////////////////////
// Dang -- 2012-056-13
// Function description:
//       tbc...
//
long ReCompXArray () {

  long k;
  ZZ x, xSum, invPowRatio;

  xSum = 0;

  // for each big set B_k
  for (k = 1; k <= SS_SIZE; k++) {
    if (k != SS_SIZE) {
      // pick a random number up to D, this number is X_(k,i_k)
      x = RandomBnd (D);
      AddMod (xSum, xSum, x, D);
    }
    else {
      // the last number X_(k,i_k) is not random
      // so that the sum of all X_(k,i_k) is exactly S (private key) mod D
      SubMod (x, S, xSum, D);
      //      x = S - xSum;

    }

    // compute X_k from X_(k,i_k)
    InvMod (invPowRatio, RATIO_POWERS(idxArray[k-1]), D);
    X_ARRAY(k) = MulMod (x, invPowRatio, D);
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-056-14
// Function description:
//       tbc...
//
long ReCompEta () {

  long j, k;
  ZZ encOne, encZero;

  // for each big set B_k
  for (k = 1; k <= SS_SIZE; k++) {

    // fill ETA matrix with encryptions of 0 or 1 (only two instances of 1)
    for (j = 1; j <= IDX_MAX; j++) {
      if (j == aArray[k-1] || j == bArray[k-1]) {
	EncryptBit (encOne, 1);
	ETA (k, j) = encOne;
      }
      else {
	EncryptBit (encZero, 0);
	ETA (k, j) = encZero;
      }
    }

  }

  return 0;
}
#endif


//////////////////////////////////////////////////////
// Dang -- 2012-05-22
// Function description:
//       tbc...
//
long CompEtaAndXArray () {

  long a, b, i, k, j;
  ZZ x, xSum, invPowRatio;
  ZZ encOne, encZero;

  idxArray.resize (SS_SIZE);

#ifdef FHE_TEST
  aArray.resize (SS_SIZE);
  bArray.resize (SS_SIZE);
#endif

  xSum = 0;

  // for each big set B_k
  for (k = 1; k <= SS_SIZE; k++) {

    do {
      // pick a random position for the first 1
      a = RandomBnd (A_MAX) + 1;

      // pick another random position for the second 1
      b = RandomBnd (IDX_MAX) + 1;

      // compute the one-dim index from a, b -> i_k = i(a,b)
      if (a < b)
	GetIndex (i, a, b);
      else
	GetIndex (i, b, a);
    } while ((a == b) || (i > BS_SIZE));

    idxArray[k - 1] = i;

#ifdef FHE_TEST
    if (a < b) {
      aArray[k - 1] = a;
      bArray[k - 1] = b;
    }
    else {
      aArray[k - 1] = b;
      bArray[k - 1] = a;
    }
#endif

    if (k != SS_SIZE) {
      // pick a random number up to D, this number is X_(k,i_k)
      x = RandomBnd (D);
      AddMod (xSum, xSum, x, D);
    }
    else {
      // the last number X_(k,i_k) is not random
      // so that the sum of all X_(k,i_k) is exactly S (private key) mod D
      SubMod (x, S, xSum, D);
      //      x = S - xSum;

    }

    // compute X_k from X_(k,i_k)
    InvMod (invPowRatio, RATIO_POWERS(i), D);
    X_ARRAY(k) = MulMod (x, invPowRatio, D);

    // fill ETA matrix with encryptions of 0 or 1 (only two instances of 1)
    for (j = 1; j <= IDX_MAX; j++) {
      if (j == a || j == b) {
	EncryptBit (encOne, 1);
	ETA(k, j) = encOne;
      }
      else {
	EncryptBit (encZero, 0);
	ETA(k, j) = encZero;
      }
    }
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-07
// Function description:
//       tbc...
//
long SlavePostProcess () {

  ZZ aCipher;
  int i, k, rank, t;
  ZZ x, y;
  ZZ y_vector, z_vector;
  RR y_prec, z_prec;

  t = (1 << PREC_PARAM);

  // receive the ciphertext from the master
  MPI_ZZ_Recv (aCipher, 0, 0);

#ifdef FHE_MPI_LOG
  cout << "[slave-" << FHE_MPI_RANK << ".SlavePostProcess] Received ciphertext from master." << endl;
#endif

  do {
    // receive job number
    FHE_MPI_Recv_Command (k, rank);

#ifdef FHE_MPI_LOG
    cout << "[slave-" << FHE_MPI_RANK << ".SlavePostProcess] Received job k=" << k << " from master." << endl;
#endif

    // if there is a job to do
    if (k) {

      clear (y_vector);
      clear (z_vector);

      // compute y_vector and z_vector
      for (i = 1; i <= BS_SIZE; i++) {

	// compute X_(k,i)
	MulMod (x, X_ARRAY(k), RATIO_POWERS(i), D);

	// compute Y_(k,i)
	MulMod (y, aCipher, x, D);
	if (IsOdd (y))
	  SetBit (y_vector, i);
	conv (y_prec, y);

	// compute Z_(k,i)
	DivPrec (z_prec, y_prec, D_PREC, PREC_PARAM);
	z_vector |= (TruncToZZ (z_prec * t) << ((PREC_PARAM + 1) * i));

      }  // end for (i=1; i <= BS_SIZE; i++)

      // slave sends a command to the master signaling that the result vectors will follow shortly
      FHE_MPI_Send_Command (k, 0);

#ifdef FHE_MPI_LOG
      cout << "[slave-" << FHE_MPI_RANK << ".SlavePostProcess] k=" << k << " has been processed. Send signal to master." << endl;
#endif

      // slave sends y_vector to master
      MPI_ZZ_Send (y_vector, 0, k);

#ifdef FHE_MPI_LOG
      cout << "[slave-" << FHE_MPI_RANK << ".SlavePostProcess] k=" << k << " has been processed. Send y_vector to master." << endl;
#endif

      // slave sends z_vector to master
      MPI_ZZ_Send (z_vector, 0, k);

#ifdef FHE_MPI_LOG
      cout << "[slave-" << FHE_MPI_RANK << ".SlavePostProcess] k=" << k << " has been processed. Send z_vector to master." << endl;
#endif

    }  // end if (k)

  } while (k);  // while there is still job to do

#ifdef FHE_MPI_LOG
  cout << "[slave-" << FHE_MPI_RANK << ".SlavePostProcess] Ending post process." << endl;
#endif

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-06
// Function description:
//       tbc...
//
long PostProcess (vec_ZZ& y_array, vec_ZZ& z_array, const ZZ& aCipher) {

  int i, k, k_done, completed_count;
  ZZ y_vector, z_vector;

  y_array.SetLength (SS_SIZE);
  z_array.SetLength (SS_SIZE);
  clear (y_array);
  clear (z_array);

  clear (y_vector);
  clear (z_vector);

  // master sends command POST_PROCESS, along with the ciphertext
  k = 1;
  for (i = 1; i < FHE_MPI_SIZE; i++) {

    // command slaves to enter PostProcess procedure
    FHE_MPI_Send_Command (FHE_MPI_COMMAND_POSTPROCESS, i);

    // send ciphertext to slaves
    MPI_ZZ_Send (aCipher, i, 0);

#ifdef FHE_MPI_LOG
    cout << "[PostProcess] Send ciphertext to " << i << endl;
#endif

    // send job k (possibly zero) to slave
    FHE_MPI_Send_Command (k, i);

#ifdef FHE_MPI_LOG
    cout << "[PostProcess] Send k=" << k << " to " << i << endl;
#endif

    // next job, otherwise k=0
    if ((k < SS_SIZE) && (k != 0))
      k++;
    else
      k = 0;

  }

  completed_count = 0;
  while (completed_count < SS_SIZE) {

    // master receives signal from slave indicating that a row k has been done
    FHE_MPI_Recv_Command (k_done, i);

#ifdef FHE_MPI_LOG
    cout << "[PostProcess] Receive signal k=" << k_done << " from " << i << endl;
#endif

    // master receives the y_vector
    MPI_ZZ_Recv (y_vector, i, k_done);
    y_array(k_done) = y_vector;

#ifdef FHE_MPI_LOG
    cout << "[PostProcess] k=" << k_done << ", received y_vector from " << i << endl;
#endif

    // master receives the z_vector
    MPI_ZZ_Recv (z_vector, i, k_done);
    z_array(k_done) = z_vector;

#ifdef FHE_MPI_LOG
    cout << "[PostProcess] k=" << k_done << ", received z_vector from " << i << endl;
#endif

    // send job k (possibly zero) to slave
    FHE_MPI_Send_Command (k, i);

#ifdef FHE_MPI_LOG
    cout << "[PostProcess] Send k=" << k << " to " << i << endl;
#endif

    // next job, otherwise k=0
    if ((k < SS_SIZE) && (k != 0))
      k++;
    else
      k = 0;

    completed_count++;
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-03
// Function description:
//       tbc...
//
long CompCumulVector (vec_ZZ& aCumulVector, const vec_ZZ& aBaseVector,
		      const long& rank) {

  long i, k;

  k = aCumulVector.length ();

  if (k > rank) {
    MulMod (aCumulVector(k), aCumulVector(k), aBaseVector(k - rank), D);
    for (i = k - 1; i > rank; i--) {
      AddMod (aCumulVector(i),
	      aCumulVector(i + 1), MulMod (aCumulVector(i), aBaseVector(i - rank), D),
	      D);
    }
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-05-30
// Function description:
//       tbc...
//
long HomGradSchlAddAndRound (ZZ& result, const mat_ZZ& Qk) {

  vec_vec_ZZ carryBits;
  vec_ZZ aBaseVector, aCumulVector;
  ZZ carryBit;

  long i, j, k, m;

  result = 0;
  carryBits.SetLength (PREC_PARAM);

  // calculate the carry bit generated by adding column i
  for (i = PREC_PARAM + 1; i > 1; i--) {

#ifdef FHE_LOG
    cout << "[HomGradSchlAddAndRound] Compute carry bit from " << i << endl;
#endif

    aBaseVector.SetLength (0);
    clear (aBaseVector);
    clear (aCumulVector);

    // first, build the base vector
    for (j = 1; j <= SS_SIZE; j++) {
      append (aBaseVector, Qk(j, i));
    }

    // append the carry bits, if any
    if (i <= PREC_PARAM) {
      append (aBaseVector, carryBits(i));
    }

    // then build the initial cumulating vector
    aCumulVector.SetLength (aBaseVector.length ());
    aCumulVector(aCumulVector.length ()) = aBaseVector(aBaseVector.length ());
    for (j = aCumulVector.length () - 1; j > 1; j--) {
      AddMod (aCumulVector(j), aBaseVector(j), aCumulVector(j + 1), D);
    }

    // the sum of the 2nd base vector will be part of the result
    if (i == 2)
      AddMod (result, aBaseVector(1), aCumulVector(2), D);

    // compute carry from column i to column j
    k = 1;
    for (j = i - 1; (j >= 1) && (aBaseVector.length () >= (m = (1 << (i - j)))); j--) {
      while (k < m) {
	CompCumulVector (aCumulVector, aBaseVector, k);
	k++;
      }

#ifdef FHE_LOG
      cout << "[HomGradSchlAddAndRound] Carry bit from " << i << " to " << j << endl;
#endif

      append (carryBits(j), aCumulVector(m));
    }
  }

  // add the sum of the 1st base vector (and the carry bits) to the result
  for (j = 1; j <= SS_SIZE; j++) {
    AddMod (result, result, Qk(j, 1), D);
  }
  for (j = carryBits(1).length (); j >= 1; j--) {
    AddMod (result, result, carryBits(1)(j), D);
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-09
// Function description:
//       tbc...
//
long SlaveRecrypt () {

  ZZ aCipher, y_vector, z_vector, linPart;
  vec_ZZ qk_vector;
  int k, rank;

  ZZ linPartA, linPartB;
  vec_ZZ roundPart;
  long a, b, i, j;

  qk_vector.SetLength (PREC_PARAM + 1);
  roundPart.SetLength (PREC_PARAM + 1);

  // receive the ciphertext from the master
  MPI_ZZ_Recv (aCipher, 0, 0);

#ifdef FHE_MPI_LOG
  cout << "[slave-" << FHE_MPI_RANK << ".SlaveRecrypt] Received ciphertext from master." << endl;
#endif

  do {
    // receive job number
    FHE_MPI_Recv_Command (k, rank);

#ifdef FHE_MPI_LOG
    cout << "[slave-" << FHE_MPI_RANK << ".SlaveRecrypt] Received job k=" << k << " from master." << endl;
#endif

    // if there is a job to do
    if (k) {

      // receive y_vector
      MPI_ZZ_Recv (y_vector, 0, k);

#ifdef FHE_MPI_LOG
      cout << "[slave-" << FHE_MPI_RANK << ".SlaveRecrypt] Received y_vector from master." << endl;
#endif

      // receive z_vector
      MPI_ZZ_Recv (z_vector, 0, k);

#ifdef FHE_MPI_LOG
      cout << "[slave-" << FHE_MPI_RANK << ".SlaveRecrypt] Received z_vector from master." << endl;
#endif

      clear (qk_vector);
      clear (roundPart);
      linPart = 0;
      linPartA = 0;
      linPartB = 0;

      i = 1;
      a = 1;
      b = 2;

      z_vector >>= (PREC_PARAM + 1);
      while (i <= BS_SIZE) {

	// compute the linear part
	if (bit (y_vector, i))
	  AddMod (linPartB, linPartB, ETA(k, b), D);

	// compute the rounding part
	for (j = PREC_PARAM + 1; j >= 1; j--) {
	  if (IsOdd (z_vector)) {
	    AddMod (roundPart(j), roundPart(j), ETA(k, b), D);
	  }
	  z_vector >>= 1;
	}

	i++;
	b++;
	if ((b > IDX_MAX) || (i > BS_SIZE)) {

	  // compute the linear part
	  AddMod (linPartA, linPartA, MulMod (ETA(k, a), linPartB, D), D);

	  // compute the rounding part
	  for (j = PREC_PARAM + 1; j >= 1; j--) {
	    AddMod (qk_vector(j), qk_vector(j), MulMod (ETA(k, a), roundPart(j), D), D);
	  }

	  a++;
	  b = a + 1;

	  linPartB = 0;
	  clear (roundPart);
	}
      }  // end while (i <= BS_SIZE)

      // compute the linear part
      AddMod (linPart, linPart, linPartA, D);

      // slave sends a command to the master signaling that the result vectors will follow shortly
      FHE_MPI_Send_Command (k, 0);

#ifdef FHE_MPI_LOG
      cout << "[slave-" << FHE_MPI_RANK << ".SlaveRecrypt] k=" << k << " has been processed. Send signal to master." << endl;
#endif

      // slave sends qk_vector to master
      MPI_VEC_ZZ_Send (qk_vector, 0, k);

#ifdef FHE_MPI_LOG
      cout << "[slave-" << FHE_MPI_RANK << ".SlaveRecrypt] k=" << k << " has been processed. Send qk_vector to master." << endl;
#endif

      // slave sends linPart to master
      MPI_ZZ_Send (linPart, 0, k);

#ifdef FHE_MPI_LOG
      cout << "[slave-" << FHE_MPI_RANK << ".SlaveRecrypt] k=" << k << " has been processed. Send linPart to master." << endl;
#endif

    }  // end if (k)

  } while (k);  // while there is still job to do

#ifdef FHE_MPI_LOG
  cout << "[slave-" << FHE_MPI_RANK << ".SlaveRecrypt] Ending post process." << endl;
#endif

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-05-24
// Function description:
//       tbc...
//
long RecryptBit (ZZ& result,
		 const vec_ZZ& y_array, const vec_ZZ& z_array,
		 const ZZ& aCipher) {

  mat_ZZ Qk;
  ZZ roundPart, linPart, linPartK;

  int i, k, k_done, completed_count;
  vec_ZZ qk_vector;

  if (y_array.length () != SS_SIZE ||
      z_array.length () != SS_SIZE) {

    cout << "[RecryptBit] y_array length=(" << y_array.length () << "), z_array length=(" << z_array.length () << ")" << endl;
    return -1;
  }

  Qk.SetDims (SS_SIZE, PREC_PARAM + 1);
  clear (Qk);
  clear (linPart);
  clear (roundPart);

  // master sends command RECRYPT
  // along with the ciphertext, y_vector and z_vector
  k = 1;
  for (i = 1; i < FHE_MPI_SIZE; i++) {

    // command slaves to enter PostProcess procedure
    FHE_MPI_Send_Command (FHE_MPI_COMMAND_RECRYPT, i);

    // send ciphertext to slaves
    MPI_ZZ_Send (aCipher, i, 0);

#ifdef FHE_MPI_LOG
    cout << "[Recrypt] Send ciphertext to " << i << endl;
#endif

    // send job k (possibly zero) to slave
    FHE_MPI_Send_Command (k, i);

#ifdef FHE_MPI_LOG
    cout << "[Recrypt] Send k=" << k << " to " << i << endl;
#endif

    if (k) {
      // send y_vector to slave
      MPI_ZZ_Send (y_array(k), i, k);

#ifdef FHE_MPI_LOG
      cout << "[Recrypt] k=" << k << ", send y_vector to " << i << endl;
#endif

      // send z_vector to slave
      MPI_ZZ_Send (z_array(k), i, k);

#ifdef FHE_MPI_LOG
      cout << "[Recrypt] k=" << k << ", send z_vector to " << i << endl;
#endif
    } // end if (k)

    // next job, otherwise k=0
    if ((k < SS_SIZE) && (k != 0))
      k++;
    else
      k = 0;

  } // end for (i=1; i <= FHE_MPI_SIZE; i++)

  completed_count = 0;
  while (completed_count < SS_SIZE) {

    // master receives signal from slave indicating that a row k has been done
    FHE_MPI_Recv_Command (k_done, i);

#ifdef FHE_MPI_LOG
    cout << "[Recrypt] Receive signal k=" << k_done << " from " << i << endl;
#endif

    // master receives the qk_vector
    MPI_VEC_ZZ_Recv (qk_vector, i, k_done);
    Qk(k_done) = qk_vector;

#ifdef FHE_MPI_LOG
    cout << "[Recrypt] k=" << k_done << ", received qk_vector from " << i << endl;
#endif

    // master receives the linPartK
    MPI_ZZ_Recv (linPartK, i, k_done);
    AddMod (linPart, linPart, linPartK, D);

#ifdef FHE_MPI_LOG
    cout << "[Recrypt] k=" << k_done << ", received linPart from " << i << endl;
#endif

    // send job k (possibly zero) to slave
    FHE_MPI_Send_Command (k, i);

#ifdef FHE_MPI_LOG
    cout << "[Recrypt] Send k=" << k << " to " << i << endl;
#endif

    if (k) {
      // send y_vector to slave
      MPI_ZZ_Send (y_array(k), i, k);

#ifdef FHE_MPI_LOG
      cout << "[Recrypt] k=" << k << ", send y_vector to " << i << endl;
#endif

      // send z_vector to slave
      MPI_ZZ_Send (z_array(k), i, k);

#ifdef FHE_MPI_LOG
      cout << "[Recrypt] k=" << k << ", send z_vector to " << i << endl;
#endif
    } // end if (k)

    // next job, otherwise k=0
    if ((k < SS_SIZE) && (k != 0))
      k++;
    else
      k = 0;

    completed_count++;
  }

  // homomorphically grade-school adding Qk's
  HomGradSchlAddAndRound (roundPart, Qk);

  // the recrypt result, finally
  AddMod (result, roundPart, linPart, D);

  return 0;
}
