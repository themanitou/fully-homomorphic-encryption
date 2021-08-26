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

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

#include <vector>

NTL_CLIENT


namespace Fhe
{

    using namespace std;

    class FullHom
    {

        // noise vector
        struct NoiseVector
        {
            long count;
            long *elements;
        };

        // variables to validate different functions
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

        // this variable is used for validation
        // it is calculated in the CompEtaAndXArray function
        vector<long> idxArray;  // this array contain the indexes i of the chunk X_k^i
        // in the big set B_k; k = 1 .. SS_SIZE

#ifdef FHE_TEST
        vector<long> aArray, bArray;
#endif

    public:

        // domain parameters
        ZZ          D, R, S;
        mat_ZZ      ETA;
        vec_ZZ      X_ARRAY;

        FullHom();
        ~FullHom();

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
        struct NoiseVector*    NoiseGen ();
        long               FreeNoise (struct NoiseVector*);
        long               PrintNoise (struct NoiseVector*);

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
        long               SlavePostProcess ();
        long               SlaveRecrypt ();

#ifdef FHE_DEBUG_LOG
        long               PrintAllVariables ();
#endif

#ifdef FHE_TEST
        void               FHE_Partial_Initialization ();
        long               ReCompXArray ();
        long               ReCompEta ();
#endif

    };

} // namespace Fhe

#endif /* FULLHOM_H_ */
