/**
 \file main.cpp
 \author Federico Negri <federico.negri@epfl.ch>
 \date 26-02-2016
 */


#include <Epetra_ConfigDefs.h>
#include "Taurus_config.h"

#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#if defined(HAVE_LIBGP) && defined(HAVE_EIGEN)

#include <taurus/Core/Utilities/GPmodel.hpp>
#include <lifev/core/LifeV.hpp>

using namespace std;
using namespace LifeV;

int
main ( int argc, char** argv )
{
    bool verbose (false);
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        verbose = true;
    }
#else
    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm () );
    verbose = true;
#endif

    GPmodel model;

    // load training points from file
    model.LoadTraining("GPtraining.txt");

    // generate and optimize GP model
    model.Generate();

    // evaluate GP model on a testing set and write to file
    model.Evaluate( 200, "GPtesting.m");

#ifdef HAVE_MPI
    if (verbose)
    {
        std::cout << "\nMPI Finalization\n" << std::endl;
    }
    MPI_Finalize();
#endif
    return ( EXIT_SUCCESS );
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    std::cout << "\n Please install Eigen and libgp \n";

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}

#endif // (HAVE_LIBGP) && (HAVE_EIGEN)
