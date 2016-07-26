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

// Includes for this problem
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include "boundaryConditions.hpp"
#include <taurus/Core/Geometry/RBF.hpp>
#include <taurus/Core/Utilities/VITreader.hpp>
#include <taurus/CSM/Models/StructuralSolverTierI.hpp>

using namespace LifeV;

int main ( int argc, char** argv )
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

    typedef RegionMesh< LinearTetra > mesh_Type;
    typedef FESpace< mesh_Type, MapEpetra > uSpaceStd_Type;
    typedef boost::shared_ptr< uSpaceStd_Type > uSpaceStdPtr_Type;
    typedef ETFESpace< mesh_Type, MapEpetra, 3, 3 > uSpaceETA_Type;
    typedef boost::shared_ptr< uSpaceETA_Type > uSpaceETAPtr_Type;
    typedef MatrixEpetra< Real > matrix_Type;
    typedef boost::shared_ptr< MatrixEpetra< Real > > matrixPtr_Type;
    typedef VectorEpetra vector_Type;

    // DATAFILE
    GetPot command_line ( argc, argv );
    const std::string dataFileName = command_line.follow ( "data", 2, "-f", "--file" );
    GetPot dataFile ( dataFileName );

    // LETTURA PARAMETRI
    VITreader reader;
    reader.Setup(dataFile("parameters/vit_file", "Training.vit"));
    std::vector<std::vector<double>> parameters = reader.getParameters();

    // LETTURA MESH
    boost::shared_ptr< mesh_Type > fullMeshPtr ( new mesh_Type ( Comm ) );
    MeshData meshData;
    meshData.setup (dataFile, "space_discretization");
    readMesh (*fullMeshPtr, meshData);

    // PARTIZIONAMENTO MESH
    boost::shared_ptr< mesh_Type > localMeshPtr;
    MeshPartitioner< mesh_Type > meshPart;
    meshPart.doPartition ( fullMeshPtr, Comm );
    localMeshPtr = meshPart.meshPartition();
    fullMeshPtr.reset();

    boost::shared_ptr< mesh_Type > localMeshPtr_ref ( new mesh_Type ( *localMeshPtr ) );

    // BCHandler to mark dofs
    boost::shared_ptr<BCHandler> bcHandler_markingDOFs( new BCHandler (*BCh_marking ( ) ) );

    // BCHandler to mark dofs
    boost::shared_ptr<BCHandler> problemBC ( new BCHandler (*BCh_problem () ) );

    // HDF5 files manager
    bool append_snapshots = dataFile("exporter/append_snapshots", false);
    bool append_system_snapshots = dataFile("exporter/append_system_snapshots", false);

    int numExportedJac  = 0;
    int numExportedRhsI = 0;
    int numExportedRhsE = 0;
    int numExportedSol = 0;

    EpetraExt::HDF5 HDF5_exporter( *Comm );
    if (append_system_snapshots)
    {
    	std::string h5_Jac_prefix     = dataFile("jacobian_matrix/prefix", "A");
    	std::string h5_resI_prefix    = dataFile("residual_internal_force/prefix", "F");
    	std::string h5_resE_prefix    = dataFile("residual_external_force/prefix", "F");

    	HDF5_exporter.Open(dataFile("exporter/name_output_file_system", "SystemSnapshots")+".h5");
    	HDF5_exporter.Read("info", "num_"+h5_Jac_prefix+"Snapshots", numExportedJac);
    	HDF5_exporter.Read("info", "num_"+h5_resI_prefix+"Snapshots", numExportedRhsI);
    	HDF5_exporter.Read("info", "num_"+h5_resE_prefix+"Snapshots", numExportedRhsE);

    }
    else
    {
    	HDF5_exporter.Create(dataFile("exporter/name_output_file_system", "SystemSnapshots")+".h5");
    }
    HDF5_exporter.Close();

    EpetraExt::HDF5 HDF5_exporterSol( *Comm );
    if (append_snapshots)
    {
    	std::string h5_sol_prefix     = dataFile("solution/prefix", "A");
    	HDF5_exporterSol.Open(dataFile("exporter/name_output_file_solution", "Snapshots")+".h5");
    	HDF5_exporterSol.Read("info", "num_"+h5_sol_prefix+"Snapshots", numExportedSol);
    }
    else
    {
    	HDF5_exporterSol.Create(dataFile("exporter/name_output_file_solution", "Snapshots")+".h5");
    }
    HDF5_exporterSol.Close();


    for ( int i = 0; i < parameters.size() ; ++i )
    {

    	if (verbose)
    	{
    		std::cout << "\n===============================================================";
    		std::cout << "\n==========  Solving for configuration " << i+1 << " of " << parameters.size() << "  =================";
    		std::cout << "\n===============================================================\n";
    	}

    	// MESH MOTION
        std::vector<double> parametersVector = parameters[i];

        StructuralSolverTierI CSM_tierI(dataFile, Comm);
        CSM_tierI.Setup(localMeshPtr, BCh_problem(), BCh_marking() );
        CSM_tierI.getMaterialModelPtr()->updateParameters(parametersVector[0], parametersVector[1], parametersVector[2]);
        CSM_tierI.SetExportOffset(numExportedJac, numExportedRhsI, numExportedRhsE, numExportedSol);

        if (i==0)
        	CSM_tierI.Solve(parametersVector, true, i);
        else
        	CSM_tierI.Solve(parametersVector, false, i);

        CSM_tierI.GetExportOffset(numExportedJac, numExportedRhsI, numExportedRhsE, numExportedSol);

        if (i == (parameters.size()-1))
        	CSM_tierI.WriteOffset();
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
