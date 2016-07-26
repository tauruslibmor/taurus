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
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/fem/BCManage.hpp>

// bc per marcare i dofs
#include "boundaryConditions.hpp"
#include <taurus/Core/Utilities/DOF_Extractor.hpp>

#include <taurus/Core/HyperReduction/HyperReductionVector.hpp>
#include <taurus/Core/HyperReduction/HyperReductionMatrix.hpp>
#include <taurus/Core/ReducedBasis/ReductionVector.hpp>

#include <taurus/Core/HyperReduction/MeshConnectivities.hpp>
#include <taurus/Core/HyperReduction/VectorizedToMatrices.hpp>
#include <taurus/Core/HyperReduction/BuildReducedMesh.hpp>


#include <taurus/Core/Utilities/Utilities.hpp>
#include <taurus/Core/Utilities/DenseHDF5.hpp>
#include <taurus/Core/Utilities/VITreader.hpp>
#include "EpetraExt_HDF5.h"

#include <iterator>
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
	typedef boost::shared_ptr< uSpaceStd_Type > uSpaceStdPtr_Type;
	typedef ETFESpace< mesh_Type, MapEpetra, 3, 3 > uSpaceETA_Type;
	typedef ETFESpace< mesh_Type, MapEpetra, 3, 3 > uSpaceETAVectorial_Type;
	typedef boost::shared_ptr< uSpaceETA_Type > uSpaceETAPtr_Type;
	typedef boost::shared_ptr< uSpaceETAVectorial_Type > uSpaceETAVectorialPtr_Type;
	typedef FESpace<mesh_Type, MapEpetra>::function_Type function_Type;
	typedef MatrixEpetra< Real > matrix_Type;
	typedef boost::shared_ptr< MatrixEpetra< Real > > matrixPtr_Type;
	typedef VectorEpetra vector_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

	// DATAFILE
	GetPot command_line ( argc, argv );
	const std::string dataFileName = command_line.follow ( "data", 2, "-f", "--file" );
	GetPot dataFile ( dataFileName );

	// GENERAZIONE MESH
	boost::shared_ptr< mesh_Type > fullMeshPtr ( new mesh_Type ( Comm ) );
	MeshData meshData;
	meshData.setup (dataFile, "space_discretization");
	readMesh (*fullMeshPtr, meshData);

	// CREAZIONE CONNETTIVITA MESH
	uSpaceStdPtr_Type uFESpace_serial ( new uSpaceStd_Type ( fullMeshPtr, dataFile ("finite_element/degree", "P1"), 3, Comm ) );
	BuildReducedMesh RedMeshGenerator(uFESpace_serial);

	uFESpace_serial.reset();

	// PARTIZIONAMENTO MESH
	boost::shared_ptr< mesh_Type > localMeshPtr;
	MeshPartitioner< mesh_Type > meshPart;
	meshPart.doPartition ( fullMeshPtr, Comm );
	localMeshPtr = meshPart.meshPartition();
	//fullMeshPtr.reset();

	// SPAZI ELEMENTI FINITI
	uSpaceStdPtr_Type uFESpace ( new uSpaceStd_Type ( localMeshPtr, dataFile ("finite_element/degree", "P1"), 3, Comm ) );
	uSpaceETAPtr_Type ETuFESpace ( new uSpaceETA_Type ( localMeshPtr, & ( uFESpace->refFE() ), & ( uFESpace->fe().geoMap() ), Comm ) );

	// OGGETTI PER VISUALIZZAZIONE MESH RIDOTTA
	uSpaceStdPtr_Type FESpace_reducedMesh ( new uSpaceStd_Type ( fullMeshPtr, dataFile ("finite_element/degree", "P1"), 1, Comm ) );
	vectorPtr_Type reduced_mesh ( new vector_Type ( FESpace_reducedMesh->map(), Unique ) );
	reduced_mesh->zero();

	ExporterHDF5< mesh_Type > exporter_reduced_mesh ( dataFile, "exporter" );
	exporter_reduced_mesh.setMeshProcId( localMeshPtr, Comm->MyPID() );
	exporter_reduced_mesh.setPrefix( dataFile ("exporter/name_output_reduced_mesh", "ReducedMesh") );
	exporter_reduced_mesh.setPostDir( "./" );
	exporter_reduced_mesh.addVariable ( ExporterData< mesh_Type >::ScalarField, "reduced mesh", FESpace_reducedMesh, reduced_mesh, UInt ( 0 ) );

    // Offset componenti
    int offsetComponent = uFESpace->dof().numTotalDof();
    
    int numSolSnapshots;
	EpetraExt::HDF5 HDF5_importer(*Comm);

	HDF5_importer.Open( dataFile ("name_output_file_system/name_input_file_solution", "SolSnapshots")+".h5" );
	std::string h5_sol_prefix     = dataFile("solution/prefix", "sol");
	HDF5_importer.Read("info", "num_"+h5_sol_prefix+"Snapshots", numSolSnapshots);
	HDF5_importer.Close();

	// ---------------------------
	// HDF5 file managers
	// ---------------------------
	EpetraExt::HDF5 HDF5_exporter( *Comm );
	HDF5_exporter.Create( dataFile("rom_exporter/name_output_file", "ROM")+".h5" );

	DenseHDF5 HDF5dense_exporter;
	if ( Comm->MyPID() == 0 )
	{
		HDF5dense_exporter.Create( dataFile("rom_exporter/hyper_output_file", "Hyper_ROM")+".h5" );
	}

	std::string JacobianApproximation = dataFile ("hyperreduction/jacobian_approximation", "DEIM");
	bool SplitInternalForces   = dataFile ("hyperreduction/split_internal_force", false);

	// ---------------------------
	// Build FE Matrix Map rows in the internal dofs
	// ---------------------------

	boost::shared_ptr<BCHandler> bc_markingDOFs = BCh_marking ();
	bc_markingDOFs->bcUpdate( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

	boost::shared_ptr<DOF_Extractor> dofExtractor;
	dofExtractor.reset ( new DOF_Extractor ( uFESpace ) );
	dofExtractor->setup( bc_markingDOFs );

	mapPtr_Type map_rows;
	map_rows.reset ( new MapEpetra ( *dofExtractor->getMapUnmarked() ) );

	// ---------------------------
	// Build Reduced Basis for the solution
	// ---------------------------
    ReductionVector RBsol(dataFile, "solution", Comm, map_rows);
	RBsol.perform();
	RBsol.write(HDF5_exporter, HDF5dense_exporter);

	boost::shared_ptr<Epetra_FECrsMatrix> V;
	RBsol.getPODbasis(V);
	int numSolBases = RBsol.getNumBases();


	// ---------------------------
	// Hyper Reduce RHS Internal ISO Forces
	// ---------------------------

	HyperReductionVector HyperFIiso( dataFile, "residual_internal_force_iso", Comm, map_rows);

	HyperFIiso.perform();
	HyperFIiso.write(HDF5dense_exporter);
	HyperFIiso.project(V, numSolBases, HDF5_exporter);

	std::vector<int> deim_FIiso_GID = HyperFIiso.getDeimGID();

	if (JacobianApproximation=="DEIM")
	{
		HyperFIiso.DEIMjacobianLeftProjection(V, HDF5_exporter, "DEIMjacobianLeftProjection_iso");
	}

	// ---------------------------
	// Hyper Reduce RHS Internal VOL Forces
	// ---------------------------
	HyperReductionVector HyperFIvol( dataFile, "residual_internal_force_vol", Comm, map_rows);

	HyperFIvol.perform();
	HyperFIvol.write(HDF5dense_exporter);
	HyperFIvol.project(V, numSolBases, HDF5_exporter);

	std::vector<int> deim_FIvol_GID = HyperFIvol.getDeimGID();

	if (JacobianApproximation=="DEIM")
	{
		HyperFIvol.DEIMjacobianLeftProjection(V, HDF5_exporter, "DEIMjacobianLeftProjection_vol");
	}

	// ---------------------------
	// Hyper Reduce RHS External Forces
	// ---------------------------
	HyperReductionVector HyperFE( dataFile, "residual_external_force", Comm, map_rows);

	HyperFE.perform();
	HyperFE.write(HDF5dense_exporter);
	HyperFE.project(V, numSolBases, HDF5_exporter);

	std::vector<int> deim_FE_GID = HyperFE.getDeimGID();

	// ---------------------------
	// Hyper Reduce MATRIX
	// ---------------------------
	if (JacobianApproximation=="MDEIM")
	{
		HyperReductionMatrix HyperA( dataFile, "jacobian_matrix", Comm, map_rows);

		HyperA.perform();
		HyperA.write(HDF5dense_exporter);
		HyperA.project(V, numSolBases, HDF5_exporter);

		std::vector<int> deim_A_GID = HyperA.getDeimDofs();
		RedMeshGenerator.appendSelectedDofs(deim_A_GID);
	}

	// ---------------------------
	// Build Reduced Mesh
	// ---------------------------
	RedMeshGenerator.appendSelectedDofs(deim_FIiso_GID);
	RedMeshGenerator.appendSelectedDofs(deim_FIvol_GID);
	RedMeshGenerator.appendSelectedDofs(deim_FE_GID);
	RedMeshGenerator.DofsToNodes_CSM();
	RedMeshGenerator.build(FESpace_reducedMesh, reduced_mesh, HDF5dense_exporter);

	exporter_reduced_mesh.postProcess(0.0);
	exporter_reduced_mesh.closeFile();
     
	// close HDF5 file manager
	HDF5_exporter.Close();
	if ( Comm->MyPID() == 0 )
	{
		HDF5dense_exporter.Close();
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
