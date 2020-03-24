/*! \file mpivar.h
 *  \brief this file declares variables used for mpi compilation.
 */

#ifndef MPIVAR_H
#define MPIVAR_H

#include <mpi.h>

///\name Include for NBodyFramework library.
//@{
///nbody code
#include <NBody.h>
///Math code
#include <NBodyMath.h>
///Binary KD-Tree code
#include <KDTree.h>
//@}

using namespace std;
using namespace Math;
using namespace NBody;


/// \name definitions for MPI types
//@{
#ifdef LONGINT
#define MPI_Int_t MPI_LONG_LONG_INT
#define MPI_UInt_t MPI_UNSIGNED_LONG_LONG
#else
#define MPI_Int_t MPI_INT
#define MPI_UInt_t MPI_UNSIGNED
#endif
#ifdef SINGLEPRECISION
#define MPI_Real_t MPI_FLOAT
#else
#define MPI_Real_t MPI_DOUBLE
#endif

#ifdef LONGIDS
#define MPI_Id_type MPI_UNSIGNED_LONG_LONG
#else
#define MPI_Id_type MPI_UNSIGNED
#endif

#ifdef TREEFROGLONGIDS
typedef long long IDTYPE;
#define MPI_IDTYPE MPI_LONG_LONG
#elif defined(TREEFROGLONGUIDS)
typedef unsigned long long IDTYPE;
#define MPI_IDTYPE MPI_UNSIGNED_LONG_LONG
#elif defined(TREEFROGINTIDS)
typedef int IDTYPE;
#define MPI_IDTYPE MPI_INTEGER
#else
typedef int unsigned IDTYPE;
#define MPI_IDTYPE MPI_UNSIGNED_INTEGER
#endif

//@}

/// \name for mpi tasks and domain construction
//@{
extern int ThisTask, NProcs;
//@}
///number of snapshots read and proccessed by thread
extern int NSnap,StartSnap,EndSnap;
extern int *mpi_startsnap,*mpi_endsnap;
extern IDTYPE StartID, EndID;
extern Int_t NumIDs;
extern vector<Int_t> mpi_numids;
extern vector<IDTYPE> mpi_startid, mpi_endid;


#endif
