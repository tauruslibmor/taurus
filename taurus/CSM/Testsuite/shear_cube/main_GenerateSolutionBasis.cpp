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

#include <taurus/Core/ReducedBasis/V_POD.hpp>
#include <taurus/Core/ReducedBasis/ReductionVector.hpp>

#include <taurus/Core/Utilities/Utilities.hpp>

#include "Epetra_Vector.h"

#include "EpetraExt_HDF5.h"
#include <iterator>
#include <chrono>



#include <taurus/Core/Utilities/DenseHDF5.hpp>
#include "hdf5.h"

#include <taurus/Core/Utilities/VITreader.hpp>

using namespace LifeV;
using namespace std::chrono;

void spy ( std::string const& fileName, boost::shared_ptr<Epetra_CrsMatrix> M_epetraCrs )
{
	std::string name = "", uti = " , ";

	Int  me = M_epetraCrs->Comm().MyPID();
	std::ostringstream myStream;
	myStream << me;
	name = fileName + ".m";

	EpetraExt::RowMatrixToMatlabFile ( name.c_str(), *M_epetraCrs );
}

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

	// PARTIZIONAMENTO MESH
	boost::shared_ptr< mesh_Type > localMeshPtr;
	MeshPartitioner< mesh_Type > meshPart;
	meshPart.doPartition ( fullMeshPtr, Comm );
	localMeshPtr = meshPart.meshPartition();
	//fullMeshPtr.reset();

	// SPAZI ELEMENTI FINITI
	uSpaceStdPtr_Type uFESpace ( new uSpaceStd_Type ( localMeshPtr, dataFile ("finite_element/degree", "P1"), 3, Comm ) );
	uSpaceETAPtr_Type ETuFESpace ( new uSpaceETA_Type ( localMeshPtr, & ( uFESpace->refFE() ), & ( uFESpace->fe().geoMap() ), Comm ) );

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
	RBsol.write(HDF5_exporter);

	// close HDF5 file manager
	HDF5_exporter.Close();



#ifdef HAVE_MPI
	if (verbose)
	{
		std::cout << "\nMPI Finalization\n" << std::endl;
	}
	MPI_Finalize();
#endif
	return ( EXIT_SUCCESS );
}
