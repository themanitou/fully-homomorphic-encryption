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

#include "calculation/fullhom.h"


//////////////////////////////////////////////////////
// Dang -- 2012-05-31
// Function description:
//       tbc...
//
long validate_Eta () {

    long i1, i2, k, a, b;
    long plainText;
    ZZ cipherText;
    int valid;

    valid = 1;

    for (k = 1; k <= SS_SIZE; k++) {
        i1 = idxArray[k-1];
        cout << "[worker-" << FHE_MPI_RANK << ".validate_Eta] idxArray[" << (k-1) << "]=" << i1 << endl;
        for (a = 1; a <= A_MAX; a++) {
            for (b = a+1; b <= IDX_MAX; b++) {
                MulMod (cipherText, ETA(k, a), ETA(k, b), D);
                DecryptBit (plainText, cipherText);
                if (plainText) {
                    GetIndex (i2, a, b);
                    cout << "[worker-" << FHE_MPI_RANK << ".validate_Eta] found for idxArray[" << (k-1) << "]::" << i2 << endl;

                    if (i2 != i1)
                        valid = 0;
                }
            }
        }
    }

    if (valid)
        cout << "[worker-" << FHE_MPI_RANK << ".validate_Eta] Success." << endl;
    else
        cout << "[worker-" << FHE_MPI_RANK << ".validate_Eta] Failed." << endl;

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-05-31
// Function description:
//       tbc...
//
long validate_Xarray () {

    long i, k;
    ZZ ratio, powRatio, sumX, diff;

    ratio = power2_ZZ (LOG2_RATIO);
    powRatio = 0;
    sumX = 0;

    for (k = 1; k <= SS_SIZE; k++) {
        i = idxArray[k - 1];
        PowerMod (powRatio, ratio, i, D);
        AddMod (sumX, sumX, MulMod (X_ARRAY(k), powRatio, D), D);
    }

    SubMod (diff, S, sumX, D);
    if (IsZero (diff))
        cout << "[worker-" << FHE_MPI_RANK << ".validate_Xarray]: Success." << endl;
    else {
        cout << "[worker-" << FHE_MPI_RANK << ".validate_Xarray]: Failed." << endl;
        cout << "[worker-" << FHE_MPI_RANK << ".validate_Xarray] S = " << S << endl;
        cout << "[worker-" << FHE_MPI_RANK << ".validate_Xarray] sumX = " << sumX << endl;
        cout << "[worker-" << FHE_MPI_RANK << ".validate_Xarray] S - sumX = " << diff << endl;
    }

    return 0;

}


//////////////////////////////////////////////////////
// Dang -- 2012-05-31
// Function description:
//       tbc...
//
long validate_RecryptBit () {
    long count, i;
    long plainText, decryptText;
    ZZ cipherText, recryptText;
    vec_ZZ y_array, z_array;

    long similar, nonsimilar, error;
    double simRate, nonSimRate, errorRate;

    clock_t ppClockTics, ppClStart, ppClStop;
    clock_t rpClockTics, rpClStart, rpClStop;
    double ppTime, rpTime;
    double totalTime;

    count = 4;
    similar = 0;
    nonsimilar = 0;
    error = 0;
    simRate = 0.0;
    nonSimRate = 0.0;
    errorRate = 0.0;

    ppClockTics = 0.0;
    rpClockTics = 0.0;

    for (i = 0; i < count; i++) {
        cout << "[validate_RecryptBit] #" << i;

        plainText = RandomBnd (2);
        decryptText = 0;

        EncryptBit (cipherText, plainText);

        ppClStart = clock ();
        PostProcess (y_array, z_array, cipherText);
        ppClStop = clock ();
        ppClockTics = (ppClStop - ppClStart);

        rpClStart = clock ();
        RecryptBit (recryptText, y_array, z_array, cipherText);
        rpClStop = clock ();
        rpClockTics = (rpClStop - rpClStart);

        DecryptBit (decryptText, recryptText);
        cout << ": [plaintext,recrypt,decrypted]=[" << plainText << "," << IsOdd (recryptText) << "," << decryptText << "]." << endl;

        if (plainText ^ IsOdd (recryptText))
            nonsimilar++;
        else
            similar++;

        if (plainText ^ decryptText)
            error++;

#ifdef FHE_DEBUG_LOG
        if ((error == 1) && (plainText ^ decryptText)) {
            PrintAllVariables ();
        }
#endif

        simRate = (double)similar/(double)(i+1);
        nonSimRate = (double)nonsimilar/(double)(i+1);
        errorRate = (double)error/(double)(i+1);

        cout << "[validate_RecryptBit] Similar=" << similar << "(" << simRate << ")";
        cout << ", Non-similar=" << nonsimilar << "(" << nonSimRate << ")";
        cout << ", Error=" << error << "(" << errorRate << ")" << endl;

        ppTime = (double)ppClockTics / (double)CLOCKS_PER_SEC;
        rpTime = (double)rpClockTics / (double)CLOCKS_PER_SEC;
        totalTime = ppTime + rpTime;

        cout << "[validate_RecryptBit] Post-process time=" << ppTime << "s (" << (ppTime/totalTime) << "). Recrypt time=" << rpTime << "s (" << (rpTime/totalTime) << "). Total time=" << totalTime << "s." << endl;

    }

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-03-06
// Function description:
//       tbc...
//
long validate_MulHalfByte () {

    long count, i;
    long plaintext1, plaintext2, plaintext3;
    long n1, n2, n3, product;
    vec_ZZ ciphertext1, ciphertext2, ciphertext3;

    long error;
    double errorRate;

    count = 200;
    error = 0;
    errorRate = 0.0;

    for (i = 0; i < count; i++) {

        cout << "[validate_MulHalfByte] #" << i;

        plaintext1 = RandomBnd (16);
        plaintext2 = RandomBnd (16);

        n1 = plaintext1;
        n2 = plaintext2;

        product = n1 * n2;

        EncryptByte (ciphertext1, plaintext1);
        EncryptByte (ciphertext2, plaintext2);
        ciphertext1[4] = ciphertext1[5] = ciphertext1[5] = ciphertext1[7] = 0;
        ciphertext2[4] = ciphertext2[5] = ciphertext2[5] = ciphertext2[7] = 0;

        MulHalfByte (ciphertext3, ciphertext1, ciphertext2);

        DecryptByte (plaintext3, ciphertext3);
        n3 = plaintext3;

        cout << ": [n1,n2,*]=[" << n1 << "," << n2 << "," << n3 << "].";

        if (n3 != product) {
            cout << "***";
            error++;
        }

        cout << endl;
    }

    errorRate = (double)error/(double)count;
    cout << "[validate_MulHalfByte] Error=" << error << "(" << errorRate << ")." << endl;

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-24
// Function description:
//       tbc...
//
long validate_MulByte () {

    long count, i;
    long plaintext1, plaintext2, plaintext3;
    long n1, n2, n3, product;
    vec_ZZ ciphertext1, ciphertext2, ciphertext3;

    long error;
    double errorRate;

    count = 100;
    error = 0;
    errorRate = 0.0;

    for (i = 0; i < count; i++) {

        cout << "[validate_MulByte] #" << i;

        plaintext1 = RandomBnd (16);
        plaintext2 = RandomBnd (16);

        n1 = plaintext1;
        n2 = plaintext2;

        product = n1 * n2;

        EncryptByte (ciphertext1, plaintext1);
        EncryptByte (ciphertext2, plaintext2);
        ciphertext1[4] = ciphertext1[5] = ciphertext1[5] = ciphertext1[7] = 0;
        ciphertext2[4] = ciphertext2[5] = ciphertext2[5] = ciphertext2[7] = 0;

        MulByte (ciphertext3, ciphertext1, ciphertext2);

        DecryptWord (plaintext3, ciphertext3);
        n3 = plaintext3;

        cout << ": [n1,n2,*]=[" << n1 << "," << n2 << "," << n3 << "].";

        if (n3 != product) {
            cout << "***";
            error++;
        }

        cout << endl;
    }

    errorRate = (double)error/(double)count;
    cout << "[validate_MulByte] Error=" << error << "(" << errorRate << ")." << endl;

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-24
// Function description:
//       tbc...
//
long validate_MaxMinSignedByte () {

    long count, i;
    long plaintext1, plaintext2, plaintext3, plaintext4;
    long n1, n2, n3, n4, bMax, bMin;
    vec_ZZ ciphertext1, ciphertext2, ciphertext3, ciphertext4;

    ZZ carry1, carry2;

    long errormax, errormin;
    double errorMaxRate, errorMinRate;

    count = 200;
    errormax = 0;
    errormin = 0;
    errorMaxRate = 0.0;
    errorMinRate = 0.0;

    for (i = 0; i < count; i++) {

        cout << "[validate_MaxMinSignedByte] #" << i;

        plaintext1 = RandomBnd (256);
        plaintext2 = RandomBnd (256);

        if (plaintext1 > 127)
            n1 = plaintext1 - 256;
        else
            n1 = plaintext1;

        if (plaintext2 > 127)
            n2 = plaintext2 - 256;
        else
            n2 = plaintext2;

        bMax = max(n1,n2);
        bMin = min(n1,n2);

        EncryptByte (ciphertext1, plaintext1);
        EncryptByte (ciphertext2, plaintext2);

        MaxMinSignedByte (ciphertext3, ciphertext4, ciphertext1, ciphertext2);

        DecryptByte (plaintext3, ciphertext3);
        if (plaintext3 > 127)
            n3 = plaintext3 - 256;
        else
            n3 = plaintext3;

        DecryptByte (plaintext4, ciphertext4);
        if (plaintext4 > 127)
            n4 = plaintext4 - 256;
        else
            n4 = plaintext4;

        cout << ": [n1,n2,max,min]=[" << n1 << "," << n2 << "," << n3 << "," << n4 << "]." << endl;

        if (n3 != bMax)
            errormax++;

        if (n4 != bMin)
            errormin++;
    }

    errorMaxRate = (double)errormax/(double)count;
    errorMinRate = (double)errormin/(double)count;
    cout << "[validate_MaxMinSignedByte] [ErrorMax, ErrorMin]=[" << errormax << "(" << errorMaxRate << ")," << errormin << "(" << errorMinRate << ")]." << endl;

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long validate_NegByte () {

    long count, i, j, b, nb;
    vec_ZZ x;
    ZZ c;

    long error;

    count = 10;

    error = 0;

    for (i = 0; i < count; i++) {
        b = RandomBnd (256);

        EncryptByte (x, b);
        NegByte (x, x);
        DecryptByte (nb, x);

        cout << "[validate_NegByte] [b,neg_b]=[" << b << "," << nb-256 << "],[";
        for (j = 7; j >= 0; j--) {
            cout << ((b >> j) % 2);
        }
        cout << ",";
        for (j = 7; j >= 0; j--) {
            cout << ((nb >> j) % 2);
        }
        cout << "]." << endl;

        if (b + nb != 256)
            ++error;

    }

    double errorRate = (double)error/(double)count;
    cout << "[validate_NegByte] Error=" << error << "(" << errorRate << ")." << endl;

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long validate_FlipBit () {

    long count, i, j, k;
    ZZ x;

    long error;

    count = 1000;

    error = 0;

    for (i = 0; i < count; i++) {
        j = RandomBnd (2);

        EncryptBit (x, j);
        FlipBit (x, x);
        DecryptBit (k, x);

        cout << "[validate_FlipBit] [b,b_flip]=[" << j << "," << k << "]." << endl;

        if (j == k)
            error++;
    }

    cout << "[validate_FlipBit] Error=" << error << endl;

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-24
// Function description:
//       tbc...
//
long validate_XorBit () {

    long count, i, b1, b2, b3, bXor;
    ZZ x1, x2, x3;

    long error;
    double errorRate;

    count = 1000;

    error = 0;
    errorRate = 0.0;

    for (i = 0; i < count; i++) {
        b1 = RandomBnd (2);
        b2 = RandomBnd (2);
        bXor = b1 ^ b2;

        EncryptBit (x1, b1);
        EncryptBit (x2, b2);
        XorBit (x3, x1, x2);
        DecryptBit (b3, x3);

        cout << "[validate_XorBit] [b1,b2,xor]=[" << b1 << "," << b2 << "," << b3 << "]." << endl;

        if (bXor != b3)
            error++;
    }

    errorRate = (double)error/(double)count;
    cout << "[validate_XorBit] Error=" << error << "(" << errorRate << ")." << endl;

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long validate_AddByte_SubByte () {

    long count, i;
    long plaintext1, plaintext2, plaintext3, plaintext4;
    long n1, n2, n3, n4, sum, sub;
    vec_ZZ ciphertext1, ciphertext2, ciphertext3, ciphertext4;

    ZZ carry1, carry2;
    //  long carrybit1, carrybit2;

    long errorsum, errorsub;
    double errorSumRate, errorSubRate;

    count = 100;
    errorsum = 0;
    errorsub = 0;
    errorSumRate = 0.0;
    errorSubRate = 0.0;

    for (i = 0; i < count; i++) {

        cout << "[validate_AddByte_SubByte] #" << i;

        plaintext1 = RandomBnd (256);
        plaintext2 = RandomBnd (256);

        if (plaintext1 > 127)
            n1 = plaintext1 - 256;
        else
            n1 = plaintext1;

        if (plaintext2 > 127)
            n2 = plaintext2 - 256;
        else
            n2 = plaintext2;

        sum = n1 + n2;
        if (sum > 127)
            sum -= 256;
        else if (sum < -128)
            sum += 256;

        sub = n1 - n2;
        if (sub > 127)
            sub -= 256;
        else if (sub < -128)
            sub += 256;

        EncryptByte (ciphertext1, plaintext1);
        EncryptByte (ciphertext2, plaintext2);

        AddByte (ciphertext3, carry1, ciphertext1, ciphertext2);

        DecryptByte (plaintext3, ciphertext3);
        if (plaintext3 > 127)
            n3 = plaintext3 - 256;
        else
            n3 = plaintext3;
        /*
          DecryptBit (carrybit1, carry1);
          plaintext3 |= (carrybit1 << 8);
        */

        /*
          NegByte (ciphertext2, ciphertext2);
          AddByte (ciphertext4, carry2, ciphertext1, ciphertext2);
        */
        SubByte (ciphertext4, carry2, ciphertext1, ciphertext2);

        DecryptByte (plaintext4, ciphertext4);
        if (plaintext4 > 127)
            n4 = plaintext4 - 256;
        else
            n4 = plaintext4;
        /*
          DecryptBit (carrybit2, carry2);
          plaintext3 |= (carrybit1 << 8);
        */

        cout << ": [n1,n2,+,-]=[" << n1 << "," << n2 << "," << n3 << "," << n4 << "]." << endl;

        if (n3 != sum)
            errorsum++;

        if (n4 != sub)
            errorsub++;
    }

    errorSumRate = (double)errorsum/(double)count;
    errorSubRate = (double)errorsub/(double)count;
    cout << "[validate_AddByte_SubByte] [ErrorSum, ErrorSub]=[" << errorsum << "(" << errorSumRate << ")," << errorsub << "(" << errorSubRate << ")]." << endl;

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long validate_EncryptByte_DecryptByte () {

    long count, i;
    long plaintext, decryptedtext;
    vec_ZZ ciphertext;

    long error;
    double errorRate;

    count = 1000;
    error = 0;
    errorRate = 0.0;

    for (i = 0; i < count; i++) {
        cout << "[validate_EncryptByte_DecryptByte] #" << i;
        plaintext = RandomBnd (256);
        decryptedtext = 0;
        EncryptByte (ciphertext, plaintext);
        DecryptByte (decryptedtext, ciphertext);
        cout << ": [plaintext,decrypted]=[" << plaintext << "," << decryptedtext << "]." << endl;

        if (plaintext != decryptedtext)
            error++;
    }

    errorRate = (double)error/(double)count;
    cout << "[validate_EncryptByte_DecryptByte] Error=" << error << "(" << errorRate << ")" << endl;

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long validate_EncryptBit_DecryptBit () {
    long count, i;
    long plaintext, decryptedtext;
    ZZ ciphertext;

    long similar, nonsimilar, error;
    double simRate, nonSimRate, errorRate;

    count = 8000;
    similar = 0;
    nonsimilar = 0;
    error = 0;
    simRate = 0.0;
    nonSimRate = 0.0;
    errorRate = 0.0;

    for (i = 0; i < count; i++) {
        cout << "[validate_EncryptBit_DecryptBit] #" << i;
        plaintext = RandomBnd (2);
        decryptedtext = 0;
        EncryptBit (ciphertext, plaintext);
        DecryptBit (decryptedtext, ciphertext);
        cout << ": [plaintext,ciphertext,decrypted]=[" << plaintext << "," << IsOdd (ciphertext) << "," << decryptedtext << "]." << endl;

        if (plaintext ^ IsOdd (ciphertext))
            nonsimilar++;
        else
            similar++;

        if (plaintext ^ decryptedtext)
            error++;

    }

    simRate = (double)similar/(double)count;
    nonSimRate = (double)nonsimilar/(double)count;
    errorRate = (double)error/(double)count;

    cout << "[validate_EncryptBit_DecryptBit] Similar=" << similar << "(" << simRate << ")";
    cout << ", Non-similar=" << nonsimilar << "(" << nonSimRate << ")";
    cout << ", Error=" << error << "(" << errorRate << ")" << endl;

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-23
// Function description:
//       tbc...
//
long validate_NoiseGen () {

    long count, i;
    noise_vector_pt pNoise;

    count = 1000;

    for (i = 0; i < count; i++) {
        cout << "[validate_NoiseGen] Noise #" << i << endl;
        pNoise = NoiseGen ();
        PrintNoise (pNoise);
        FreeNoise (pNoise);
    }

    return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-02-01
// Function description:
//       tbc...
//
long validate_KeyGen_NaifKeyGen () {

    long count, kgSuccess, kgFailure, nkgSuccess, nkgFailure, i;
    double kgSuccessRate, kgFailureRate, nkgSuccessRate, nkgFailureRate;
    long errCount;
    double errRate;
    long kgRes, nkgRes;
    ZZ kgD, kgR, kgS, nkgD, nkgR, nkgS;
    ZZX kgV, nkgV;

    clock_t kgClockTics, nkgClockTics, clStart, clStop;
    double kgMs, nkgMs, msRatio;

    count = 100;
    kgFailure = 0;
    kgSuccess = 0;
    nkgFailure = 0;
    nkgSuccess = 0;
    errCount = 0;

    kgClockTics = 0;
    nkgClockTics = 0;

    // Generate public and private keys
    for (i = 1; i <= count; i++) {

        cout << "[validate_KeyGen_NaifKeyGen] Generate keys #" << i << ": \n";

        kgV = RandPolyGen ();
        nkgV = kgV;
        D = 0;
        R = 0;
        S = 0;

        clStart = clock ();
        kgRes = PrimaryKeyGen (kgV);
        clStop = clock ();
        kgClockTics += clStop - clStart;

        kgD = D;
        kgR = R;
        kgS = S;

        if (kgRes)
            kgFailure++;
        else
            kgSuccess++;

        switch (kgRes) {

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

        clStart = clock ();
        nkgRes = NaifPrimaryKeyGen (nkgV);
        clStop = clock ();
        nkgClockTics += clStop - clStart;

        nkgD = D;
        nkgR = R;
        nkgS = S;

        if (nkgRes)
            nkgFailure++;
        else
            nkgSuccess++;

        switch (nkgRes) {

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

        if (kgD != nkgD) {
            cout << "[validate_KeyGen_NaifKeyGen] ERROR: D differ.\n";
        }

        if (kgR != nkgR) {
            cout << "[validate_KeyGen_NaifKeyGen] ERROR: R differ.\n";
        }

        if (kgS != nkgS) {
            cout << "[validate_KeyGen_NaifKeyGen] ERROR: S differ.\n";
        }

        if (kgRes != nkgRes ||
            kgD != nkgD ||
            kgR != nkgR ||
            kgS != nkgS)
            errCount++;
    }

    kgSuccessRate = (double)kgSuccess / (double)count;
    kgFailureRate = (double)kgFailure / (double)count;
    nkgSuccessRate = (double)nkgSuccess / (double)count;
    nkgFailureRate = (double)nkgFailure / (double)count;
    errRate = (double)errCount / (double)count;

    kgMs = (double)kgClockTics / (double)CLOCKS_PER_SEC / (double)count;
    nkgMs = (double)nkgClockTics / (double)CLOCKS_PER_SEC / (double)count;
    msRatio = (double)nkgClockTics / (double)kgClockTics;

    cout << "[validate_KeyGen_NaifKeyGen] Success: KeyGen()=" << kgSuccess << " (" << kgSuccessRate << "), NaifKeyGen()="<< nkgSuccess << " (" << nkgSuccessRate <<").\n";
    cout << "[validate_KeyGen_NaifKeyGen] Failure: KeyGen()=" << kgFailure << " (" << kgFailureRate << "), NaifKeyGen()="<< nkgFailure << " (" << nkgFailureRate <<").\n";
    cout << "[validate_KeyGen_NaifKeyGen] Errors=" << errCount << " (" << errRate << ").\n";
    cout << "[validate_KeyGen_NaifKeyGen] Time: NaifKeyGen=" << nkgMs << "(s), KeyGen=" << kgMs << "(s), NaifKeyGen/KeyGen=" << msRatio << ".\n";

    return errCount;
}


//////////////////////////////////////////////////////
// Dang -- 2012-01-31
// Function description:
//       tbc...
//
int main (int argc, char *argv[]) {

    int choice, i;
    int count;

    Initialization (argc, argv);
    count = 0;
    choice = 0;

    do {

        // after initialization, manager will interact with user, while workers
        // will wait for manager's command
        if (FHE_MPI_RANK == 0) {
            cout << "[main] Fully Homomorphic Encryption Test:" << endl;

            cout << "[main] 1 - KeyGen & NaifKeyGen." << endl;
            cout << "[main] 2 - NoiseGen." << endl << endl;

            cout << "[main] 3 - EncryptBit & DecryptBit." << endl;
            cout << "[main] 4 - EncryptByte & DecryptByte." << endl << endl;

            cout << "[main] 5 - FlipBit." << endl;
            cout << "[main] 6 - XorBit." << endl << endl;

            cout << "[main] 7 - NegByte." << endl;
            cout << "[main] 8 - AddByte & SubByte." << endl << endl;

            cout << "[main] 9 - MaxMinSignedByte." << endl << endl;

            cout << "[main] 10 - Grade-school 4-bit Multiplication." << endl;
            cout << "[main] 11 - Flat 4-bit Multiplication." << endl << endl;

            cout << "[main] 12 - RecryptBit." << endl << endl;

            cout << "[main] 13 - X_ARRAY." << endl;
            cout << "[main] 14 - Eta." << endl << endl;

            cout << "[main] 15 - Re-initialization." << endl << endl;

            cout << "[main] 0 - Quit." << endl << endl;

            cout << "[main] Test #" << flush;
            cin >> choice;

            /*
            if (count < 200)
                if (choice == 12)
                    choice = 15;
                else
                    choice = 12;
            else
                choice = 0;
            */
            cout << count << "." << choice << endl;
            count++;

            switch (choice) {

            case 1:
                validate_KeyGen_NaifKeyGen ();
                break;

            case 2:
                validate_NoiseGen ();
                break;

            case 3:
                validate_EncryptBit_DecryptBit ();
                break;

            case 4:
                validate_EncryptByte_DecryptByte ();
                break;

            case 5:
                validate_FlipBit ();
                break;

            case 6:
                validate_XorBit ();
                break;

            case 7:
                validate_NegByte ();
                break;

            case 8:
                validate_AddByte_SubByte ();
                break;

            case 9:
                validate_MaxMinSignedByte ();
                break;

            case 10:
                validate_MulByte ();
                break;

            case 11:
                validate_MulHalfByte ();
                break;

            case 12:
                validate_RecryptBit ();
                break;

            case 13:
                // manager sends command VAL_XARRAY to workers
                for (i = 1; i < FHE_MPI_SIZE; i++)
                    FHE_MPI_Send_Command (FHE_MPI_COMMAND_VAL_XARRAY, i);
                break;

            case 14:
                // manager sends command VAL_ETA to workers
                for (i = 1; i < FHE_MPI_SIZE; i++)
                    FHE_MPI_Send_Command (FHE_MPI_COMMAND_VAL_ETA, i);
                break;

            case 15:
                // manager sends command RE-INITIALIZATION
                // upon receiving this command, workers will get ready to receive
                // new domain parameters
                for (i = 1; i < FHE_MPI_SIZE; i++)
                    FHE_MPI_Send_Command (FHE_MPI_COMMAND_INITIALIZATION, i);

                // manager generates new domain parameters and send them to workers
                FHE_Initialization ();
                break;

            case 0:
                // manager sends command EXIT
                for (i = 1; i < FHE_MPI_SIZE; i++)
                    FHE_MPI_Send_Command (FHE_MPI_COMMAND_EXIT, i);

            }

        } // end manager

        else { // worker

#ifdef FHE_MPI_LOG
            cout << "[worker-" << FHE_MPI_RANK << ".main] waiting for command..." << endl;
#endif

            FHE_MPI_Recv_Command (choice, i);

#ifdef FHE_MPI_LOG
            cout << "[worker-" << FHE_MPI_RANK << ".main] recv command " << choice << endl;
#endif

            switch (choice) {

            case FHE_MPI_COMMAND_POSTPROCESS: // post-process a ciphertext
                WorkerPostProcess ();
                break;

            case FHE_MPI_COMMAND_RECRYPT: // post-process a ciphertext
                WorkerRecrypt ();
                break;

            case FHE_MPI_COMMAND_INITIALIZATION: // re-init, receive new domain params
                FHE_Initialization ();
                break;

            case FHE_MPI_COMMAND_VAL_XARRAY: // validate XARRAY
                validate_Xarray ();
                break;

            case FHE_MPI_COMMAND_VAL_ETA: // validate XARRAY
                validate_Eta ();
                break;

            } // end switch (choice)

        } // end if (FHE_MPI_RANK == 0)

    } while (choice != FHE_MPI_COMMAND_EXIT);

    CloseDown ();

    return 0;
}
