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
#include <taurus/CSM/Models/StructuralSolverTierII.hpp>

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

#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER  start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name)  if (verbose) std::cout << "\n\nRUNTIME of " << name << ": " << \
		std::chrono::duration_cast<std::chrono::milliseconds>( \
				std::chrono::high_resolution_clock::now()-start \
		).count()/1000.0 << " s. \n\n";

    typedef RegionMesh< LinearTetra > mesh_Type;
    typedef FESpace< mesh_Type, MapEpetra > uSpaceStd_Type;
    typedef boost::shared_ptr< uSpaceStd_Type > uSpaceStdPtr_Type;
    typedef ETFESpace< mesh_Type, MapEpetra, 3, 3 > uSpaceETA_Type;
    typedef boost::shared_ptr< uSpaceETA_Type > uSpaceETAPtr_Type;
    typedef MatrixEpetra< Real > matrix_Type;
    typedef boost::shared_ptr< MatrixEpetra< Real > > matrixPtr_Type;
    typedef VectorEpetra vector_Type;

    INIT_TIMER

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

    uSpaceStdPtr_Type uFESpace ( new uSpaceStd_Type ( localMeshPtr, dataFile ("finite_element/degree", "P1"), 3, Comm ) );

    boost::shared_ptr< mesh_Type > localMeshPtr_ref ( new mesh_Type ( *localMeshPtr ) );

    // BCHandler to mark dofs
    boost::shared_ptr<BCHandler> bcHandler_markingDOFs( new BCHandler (*BCh_marking ( ) ) );

    // BCHandler to mark dofs
    boost::shared_ptr<BCHandler> problemBC ( new BCHandler (*BCh_problem () ) );

    // HDF5 file manager
    EpetraExt::HDF5 HDF5_exporter( *Comm );
    HDF5_exporter.Create(dataFile("exporter/name_output_file_system", "SystemSnapshots")+".h5");
    HDF5_exporter.Close();

    int numExportedJac = 0;
    int numExportedRhsI = 0;
    int numExportedRhsE = 0;

    /////////////////////////////////////////////////////
    // Load Matrix V
    DOF_Extractor extractor ( uFESpace ) ;
    extractor.setup( bcHandler_markingDOFs );

    EpetraExt::HDF5 HDF5_importer(*Comm);
    HDF5_importer.Open( dataFile ("importer/file_name", "ROM.h5") );

    int num_basis;
    HDF5_importer.Read("info", "num_basis", num_basis);

    std::vector<int> column_indices;
    for ( int i = 0; i < num_basis; ++i )
    	column_indices.push_back(i);

    int* pointerToDofs (0);
    if (column_indices.size() > 0)
    	pointerToDofs = &column_indices[0];

    boost::shared_ptr<MapEpetra> map_column_V ( new MapEpetra ( -1, static_cast<int> (column_indices.size() ), pointerToDofs, Comm ) );

    Epetra_Map map_col(*map_column_V->map(Unique));

    if (verbose)
    	std::cout << "\n Loading matrix V ... \n";

    START_TIMER;
    Epetra_CrsMatrix* V;
    HDF5_importer.Read ( "V", map_col, *extractor.getMapUnmarked()->map(Unique), V );
    HDF5_importer.Close();
    STOP_TIMER("loading matrix V");
    /////////////////////////////////////////////////////

    // START SOLVING FOR DIFFERENT CONFIGURATIONS


    for ( int i = 0; i < parameters.size() ; ++i )
    {

    	if (verbose)
    	{
    		std::cout << "\n===============================================================";
    		std::cout << "\n==========  Solving for configuration " << i+1 << " of " << parameters.size() << "  =================";
    		std::cout << "\n===============================================================\n";
    	}

    	// MESH MOTION
        std::vector<double> GeoParamValues = parameters[i];

        StructuralSolverTierII CSM_tierII(dataFile, Comm);
        CSM_tierII.Setup(localMeshPtr, BCh_problem(), BCh_marking(), V);
        CSM_tierII.getMaterialModelPtr()->updateParameters(GeoParamValues[0], GeoParamValues[1], GeoParamValues[2]);
        CSM_tierII.SetExportOffset(numExportedJac, numExportedRhsI, numExportedRhsE);

        if (i==0)
        	CSM_tierII.Solve(GeoParamValues, true, i);
        else
        	CSM_tierII.Solve(GeoParamValues, false, i);

        CSM_tierII.GetExportOffset(numExportedJac, numExportedRhsI, numExportedRhsE);

        if (i == (parameters.size()-1))
        	CSM_tierII.WriteOffset();

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
