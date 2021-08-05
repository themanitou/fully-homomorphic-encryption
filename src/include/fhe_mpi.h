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


#ifndef FHE_MPI_H_
#define FHE_MPI_H_

#include <mpi.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>


//////////////////////////////////////////////////////
// Dang -- 2012-06-07
//
// Definitions
//
#define FHE_MPI_MAX_BUFFER               50000
#define FHE_MPI_MAX_BIT_LENGTH           200000

#define FHE_MPI_COMMAND_EXIT             0
#define FHE_MPI_COMMAND_POSTPROCESS      1
#define FHE_MPI_COMMAND_RECRYPT          2
#define FHE_MPI_COMMAND_INITIALIZATION   3
#define FHE_MPI_COMMAND_VAL_XARRAY       4
#define FHE_MPI_COMMAND_VAL_ETA          5


////////////////////////////////////////////////////
// Dang -- 2012-06-07
//
//    Global variable declarations
//
extern int               FHE_MPI_SIZE;
extern int               FHE_MPI_RANK;
extern unsigned char*    FHE_MPI_BUFFER;


////////////////////////////////////////////////////
// Dang -- 2012-06-07
//
//    Function prototypes
//
long               MPI_ZZ_Bcast (NTL::ZZ&);
long               MPI_ZZ_Send (const NTL::ZZ&, const int&, const int&);
long               MPI_ZZ_Recv (NTL::ZZ&, const int&, const int&);

long               MPI_VEC_ZZ_Send (const NTL::vec_ZZ&, const int&, const int&);
long               MPI_VEC_ZZ_Recv (NTL::vec_ZZ&, const int&, const int&);

long               FHE_MPI_Send_Command (const int&, const int&);
long               FHE_MPI_Recv_Command (int&, int&);


#endif /* FHE_MPI_H_ */
