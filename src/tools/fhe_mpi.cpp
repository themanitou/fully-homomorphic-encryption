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


#include <mpi.h>
#include "../include/fhe_mpi.h"
#include "../include/fullhom.h"


////////////////////////////////////////////////////
// Dang -- 2012-06-07
//
//    Global variable declarations
//
int               FHE_MPI_SIZE;
int               FHE_MPI_RANK;
unsigned char*    FHE_MPI_BUFFER = (unsigned char*) malloc (FHE_MPI_MAX_BUFFER);


//////////////////////////////////////////////////////
// Dang -- 2012-06-07
// Function description:
//       tbc...
//
long MPI_ZZ_Bcast (NTL::ZZ& z) {

  int n;
  char signIndicator;

  // only the master prepare the buffer to broadcast
  if (FHE_MPI_RANK == 0) {
    // size of z
    n = NTL::NumBytes (z);

    // buffer of z
    NTL::BytesFromZZ (FHE_MPI_BUFFER, z, n);

    // sign indicator of z
    if (NTL::sign (z) == -1)
      FHE_MPI_BUFFER[n] = 1;
    else
      FHE_MPI_BUFFER[n] = 0;
  }

  // first broadcast the size of z
  MPI::COMM_WORLD.Bcast (&n, 1, MPI::INT, 0);

  // block everything until all processes have reached this points
  MPI::COMM_WORLD.Barrier ();

  // if the length of z is zero then z=0
  if (n == 0) {
    // if not master then just set z=0
    if (FHE_MPI_RANK)
      z = 0;
  }
  else {
    // then broadcast the buffer of z
    MPI::COMM_WORLD.Bcast (FHE_MPI_BUFFER, FHE_MPI_MAX_BUFFER, MPI::BYTE, 0);

    // block everything until all processes have reached this points
    MPI::COMM_WORLD.Barrier ();

    // if not master then reconstruct z from the buffer
    if (FHE_MPI_RANK) {
      signIndicator = FHE_MPI_BUFFER[n];
      NTL::ZZFromBytes (z, FHE_MPI_BUFFER, n);

      // look at the sign indicator
      if (signIndicator)
	z = -z;
    }
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-07
// Function description:
//       tbc...
//
long MPI_ZZ_Send (const NTL::ZZ& z, const int& dest, const int& tag) {

  long n;

  // look at the number of bytes (n) representing z
  // the n+1'th byte represent the sign indicator
  // careful though: if z=0 then n=0, in this case we only send 1 byte: zero
  n = NTL::NumBytes (z);
  NTL::BytesFromZZ (FHE_MPI_BUFFER, z, n);
  if (NTL::sign (z) == -1)
    FHE_MPI_BUFFER[n] = 1;
  else
    FHE_MPI_BUFFER[n] = 0;

  // send n+1 bytes
  MPI::COMM_WORLD.Send (FHE_MPI_BUFFER, n + 1, MPI::BYTE, dest, tag);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-07
// Function description:
//       tbc...
//
long MPI_ZZ_Recv (NTL::ZZ& z, const int& source, const int& tag) {

  long n;
  MPI::Status status;
  char signIndicator;

  // receive n+1 bytes
  MPI::COMM_WORLD.Recv (FHE_MPI_BUFFER, FHE_MPI_MAX_BUFFER,
			MPI::BYTE, source, tag, status);
  n = status.Get_count(MPI::BYTE) - 1;

  // if we receive only 1 byte then it is zero
  if (n == 0)
    z = 0;
  else {
    // look at the sign indicator
    signIndicator = FHE_MPI_BUFFER[n];
    NTL::ZZFromBytes (z, FHE_MPI_BUFFER, n);

    if (signIndicator)
      z = -z;
  }

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-08
// Function description:
//       tbc...
//
long MPI_VEC_ZZ_Send (const NTL::vec_ZZ& z_vector,
		      const int& dest, const int& tag) {

  long i, l;

  // send the length of the vector first
  l = z_vector.length ();
  MPI::COMM_WORLD.Send (&l, 1, MPI::LONG, dest, tag);

  // then send the elements of the vector, one by one
  for (i = 1; i <= l; i++)
    MPI_ZZ_Send (z_vector(i), dest, i);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-07
// Function description:
//       tbc...
//
long MPI_VEC_ZZ_Recv (NTL::vec_ZZ& z_vector,
		      const int& source, const int& tag) {

  long i, l;
  MPI::Status status;

  // receive the length of the vector
  MPI::COMM_WORLD.Recv (&l, 1, MPI::LONG, source, tag, status);

  // set the correct length for the vector
  z_vector.SetLength (l);

  // receive the elements of the vector
  for (i = 1; i <= l; i++)
    MPI_ZZ_Recv (z_vector(i), source, i);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-07
// Function description:
//       tbc...
//
long FHE_MPI_Send_Command (const int& aCommand, const int& rank) {

  MPI::COMM_WORLD.Send (&aCommand, 1, MPI::INT, rank, 0);

  return 0;
}


//////////////////////////////////////////////////////
// Dang -- 2012-06-07
// Function description:
//       tbc...
//
long FHE_MPI_Recv_Command (int& aCommand, int& rank) {

  MPI::Status status;

  MPI::COMM_WORLD.Recv (&aCommand, 1, MPI::INT, MPI_ANY_SOURCE, 0, status);
  rank = status.Get_source ();

  return 0;
}
