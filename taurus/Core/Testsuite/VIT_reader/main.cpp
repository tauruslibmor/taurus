/**
 \file main.cpp
 \author Davide Forti <davide.forti@epfl.ch>
 \date 10-11-2015
 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <taurus/Core/Utilities/VITreader.hpp>

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

    GetPot command_line ( argc, argv );
    const std::string dataFileName = command_line.follow ( "data", 2, "-f", "--file" );
    GetPot dataFile ( dataFileName );

    VITreader reader;
    reader.Setup(dataFile("parameters/vit_file", "Training.vit"));
    reader.showMe(Comm);

    //
    // Checking that the reading was ok
    //

    std::vector<std::vector<double>> parameters = reader.getParameters();
    if ( Comm->MyPID() == 0 )
    {
    	std::cout << "\n\nNumber of Configurations: " << parameters.size()
    			  << ", Number of Parameters: " << parameters[0].size() << std::endl;

    	for (int i = 0; i < parameters.size(); i++)
    	{
    		std::cout << "Configuration " << i+1 << ": ";
    		for (int j = 0; j < parameters[i].size(); j++)
    		{
    			std::cout << parameters[i][j] << "  ";
    		}
    		std::cout << "\n";
    	}
    }

#ifdef HAVE_MPI
    if (verbose)
    {
        std::cout << "\nMPI Finalization\n" << std::endl;
    }
    MPI_Finalize();
#endif
    return ( EXIT_SUCCESS );
}
