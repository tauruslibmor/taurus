#include <taurus/Core/ReducedBasis/BuildReducedBasis.hpp>

#include "EpetraExt_HDF5.h"
#include <taurus/Core/Utilities/Utilities.hpp>


//#include "Epetra_Vector.h"

using namespace std;

namespace LifeV
{

//=========================================================================
// Constructor //
BuildReducedBasis::BuildReducedBasis (const fespacePtr_Type& fe_space, bcPtr_Type bc, const GetPot& dataFile, const std::string& dataFile_section) :
	M_FESpace          ( fe_space ),
	M_comm             ( fe_space->map().commPtr() ),
	M_datafile         ( dataFile ),
	M_datafileSection  ( dataFile_section ),
	M_NumSnapshots	   ( 0 ),
	M_NumPODbases      ( 0 )
{
	M_bc_markingDOFs = bc;
	std::string data_field = M_datafileSection+"/prefix";
	M_prefix = M_datafile(data_field.c_str(), "soluzione");
}
//=========================================================================
void BuildReducedBasis::perform( const int NumSnapshots, const std::string solution_type)
{

	M_NumSnapshots = NumSnapshots;

	// READING SOLUTION SNAPSHOTS
	if ( M_comm->MyPID() == 0 )
	{
		std::cout << "\n Read " << M_prefix << " snapshots ...\n";

	}

	ImportSnapshots reader (M_FESpace, M_NumSnapshots, M_datafile, M_bc_markingDOFs );

	if (solution_type=="scalar")
	{
		reader.readSnapshots( M_prefix );
	}

	if (solution_type=="vectorial")
	{
		reader.readSnapshotsVectorial( M_prefix );
	}

	reader.getMatrix( M_Snapshots );

	// PERFORM POD ON SOLUTION SNAPSHOTS
	if ( M_comm->MyPID() == 0 )
	{
		std::cout << "\n POD on " << M_prefix << " snapshots ...\n";
	}

	V_POD<int> podSOL ( M_comm, M_Snapshots, M_NumSnapshots );

	std::string data_field = M_datafileSection+"/nortol";
	podSOL.performPOD_sol ( M_PODbasis, M_datafile(data_field.c_str(), 1e-2) );

	std::vector<double> M_SingularValues = podSOL.getSingularValues();

	M_NumPODbases = podSOL.getNumCols();

}
//=========================================================================
void BuildReducedBasis::write( EpetraExt::HDF5& HDF5_exporter )
{
	HDF5_exporter.Write("V", *M_PODbasis);// Save solution basis V
	HDF5_exporter.Write("info", "num_basis", M_NumPODbases);
}
//=========================================================================
void BuildReducedBasis::write( EpetraExt::HDF5& HDF5_exporter, DenseHDF5& HDF5dense_exporter )
{
	if ( M_comm->MyPID() == 0 )
	{
		std::string this_prefix = "/"+M_prefix+"/";
    	HDF5dense_exporter.WriteVectorDouble(this_prefix, "SingularValues", M_SingularValues);
 	}

	HDF5_exporter.Write("V", *M_PODbasis);// Save solution basis V
	HDF5_exporter.Write("info", "num_basis", M_NumPODbases);

}
//=========================================================================
void BuildReducedBasis::getPODbasis( boost::shared_ptr<Epetra_FECrsMatrix>& matrix )
{

	matrix = M_PODbasis;

}
//=========================================================================
} // end namespace LifeV
