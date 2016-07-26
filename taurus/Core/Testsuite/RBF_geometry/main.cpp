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
#include <lifev/core/filter/ExporterHDF5.hpp>

#include <taurus/Core/Geometry/RBF.hpp>
#include <taurus/Core/Utilities/VITreader.hpp>

#include <chrono>

using namespace LifeV;
using namespace std::chrono;

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

    typedef RegionMesh< LinearTetra > mesh_Type;
    typedef FESpace< mesh_Type, MapEpetra > uSpaceStd_Type;
    typedef boost::shared_ptr<uSpaceStd_Type> uSpaceStdPtr_Type;
    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    // DATAFILE
    GetPot command_line ( argc, argv );
    const std::string dataFileName = command_line.follow ( "data", 2, "-f", "--file" );
    GetPot dataFile ( dataFileName );

    // LETTURA PARAMETRI
    VITreader reader;
    reader.Setup(dataFile("parameters/vit_file", "Training.vit"));
    std::vector<std::vector<double>> parameters = reader.getParameters();

    // GENERAZIONE MESH
    boost::shared_ptr< mesh_Type > fullMeshPtr ( new mesh_Type ( Comm ) );
    std::string mesh_source = dataFile ("input_mesh/source", "file");

    std::cout << "mesh_source = " << mesh_source << "\n\n";

    if ( mesh_source == "file" )
    {
    	MeshData meshData;
    	meshData.setup (dataFile, "space_discretization");
    	readMesh (*fullMeshPtr, meshData);
    }
    else
    {
    	regularMesh3D ( *fullMeshPtr, 0, dataFile ("mesh/nx", 10), dataFile ("mesh/ny", 10), dataFile ("mesh/nz", 10), false, 2.4, 2.4, 2.4, -1.2, -1.2, -1.2 );
    }

    // PARTIZIONAMENTO MESH
    boost::shared_ptr< mesh_Type > localMeshPtr;
    MeshPartitioner< mesh_Type > meshPart;
    meshPart.doPartition ( fullMeshPtr, Comm );
    localMeshPtr = meshPart.meshPartition();
    fullMeshPtr.reset();

    // SPAZI ELEMENTI FINITI
    uSpaceStdPtr_Type displacementMesh_FESpace ( new uSpaceStd_Type ( localMeshPtr, "P1", 3, Comm ) ); // FESpace per spostamento sempre P1

    // Vector used for the simulations
    vectorPtr_Type displacement( new vector_Type ( displacementMesh_FESpace->map(), Unique ) );

    // Set Exporter
    ExporterHDF5< mesh_Type > exporter ( dataFile, "exporter" );
    exporter.setMeshProcId( localMeshPtr, Comm->MyPID() );
    exporter.setPrefix( dataFile ("exporter/name_output_file", "laplacien_ridotto") );
    exporter.setPostDir( "./" );
    exporter.addVariable ( ExporterData< mesh_Type >::VectorField, "displacement mesh", displacementMesh_FESpace, displacement, UInt ( 0 ) );
    exporter.postProcess( 0.0 );

    //===================================================
    //======== MESH DEFORMATION =========================

    int id_configuration = 0;
    RBF moveMeshRBF;
    moveMeshRBF.Setup("RBF_data");
    std::vector<double> GeoParamValues = parameters[id_configuration];
    double Radius[] = {0.4, 0.4, 0.001};

    high_resolution_clock::time_point tStart = high_resolution_clock::now();

    moveMeshRBF.Build("gaussian", GeoParamValues, Radius);
    //moveMeshRBF.MoveMesh(localMeshPtr);			  // muove mesh
    moveMeshRBF.MoveMesh(localMeshPtr, displacement); // muove mesh e esporta anche spostamento

    high_resolution_clock::time_point tEnd = high_resolution_clock::now();
    if ( Comm->MyPID() == 0 )
    	std::cout << "\n\nElapsed time for RBF mesh deformation: " << duration_cast<milliseconds>( tEnd - tStart ).count() << " milliseconds\n\n";

    //===================================================
    //===================================================
    
    // Export Deformed Mesh
    exporter.postProcess( 1.0 );
	exporter.closeFile();

	moveMeshRBF.ExportReferencePoints(dataFile ("exporter/name_output_file", "Solution"));
	moveMeshRBF.ExportDeformedPoints(dataFile ("exporter/name_output_file", "Solution"));

#ifdef HAVE_MPI
    if (verbose)
    {
        std::cout << "\nMPI Finalization\n" << std::endl;
    }
    MPI_Finalize();
#endif
    return ( EXIT_SUCCESS );
}
