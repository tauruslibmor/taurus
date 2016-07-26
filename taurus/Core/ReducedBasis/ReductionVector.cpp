#include <taurus/Core/ReducedBasis/ReductionVector.hpp>
#include <taurus/Core/HyperReduction/ImportRhsSnapshots.hpp>
#include <taurus/Core/ReducedBasis/V_POD.hpp>
#include "EpetraExt_HDF5.h"
#include <taurus/Core/Utilities/Utilities.hpp>


//#include "Epetra_Vector.h"

using namespace std;

namespace LifeV
{

//=========================================================================
// Constructor //
ReductionVector::ReductionVector (const GetPot& dataFile, const std::string& dataFile_section, const commPtr_Type& communicator, const mapPtr_Type& map_rows) :
	M_comm             ( communicator ),
	M_datafile         ( dataFile ),
	M_datafileSection  ( dataFile_section ),
	M_NumSnapshots	   ( 0 ),
	M_NumPODbases      ( 0 ),
	M_map_rows		   ( map_rows ),
	M_useWeightMatrix  (false)
{
	std::string data_field = M_datafileSection+"/prefix";
	M_prefix = M_datafile(data_field.c_str(), "sol");
}
//=========================================================================
void ReductionVector::setWeightMatrix ( boost::shared_ptr<Epetra_FECrsMatrix> matrix )
{
	M_weightMatrix = matrix;
	M_useWeightMatrix = true;
}
//=========================================================================
void ReductionVector::perform( )
{

	EpetraExt::HDF5 HDF5_importer(*M_comm);

	// Read Snapshots
	HDF5_importer.Open( M_datafile ("importer/name_input_file_solution", "SolSnapshots")+".h5" );
	HDF5_importer.Read("info", "num_"+M_prefix+"Snapshots", M_NumSnapshots);
	HDF5_importer.Close();

	ImportRhsSnapshots RHSreader (M_map_rows, M_comm, M_NumSnapshots, M_datafile, M_datafileSection );

	RHSreader.readSnapshots( );

	//M_Snapshots->zero();
	RHSreader.getMatrix( M_Snapshots );

	// POD
	if ( M_comm->MyPID() == 0 )
	{
		std::cout << "\n POD on " << M_datafileSection << "  snapshots ...\n";
	}

	V_POD<int> podRHS ( M_comm, M_Snapshots, M_NumSnapshots );


	std::string data_field = M_datafileSection+"/nortol";
	if (M_useWeightMatrix)
	{
		podRHS.setWeightMatrix(M_weightMatrix);
	}

	podRHS.performPOD_vector( M_PODbasis, M_datafile(data_field.c_str(), 1e-2) );

	M_SingularValues = podRHS.getSingularValues();

	M_NumPODbases  = podRHS.getNumCols();

}
//=========================================================================
void ReductionVector::write( EpetraExt::HDF5& HDF5_exporter, bool usePrefix )
{
	if (usePrefix)
	{
		HDF5_exporter.Write("V_"+M_prefix, *M_PODbasis);// Save solution basis V
		HDF5_exporter.Write("info", "num_basis_"+M_prefix, M_NumPODbases);
	}
	else
	{
		HDF5_exporter.Write("V", *M_PODbasis);// Save solution basis V
		HDF5_exporter.Write("info", "num_basis", M_NumPODbases);
	}
}
//=========================================================================
void ReductionVector::write( EpetraExt::HDF5& HDF5_exporter, DenseHDF5& HDF5dense_exporter, bool usePrefix )
{
	if ( M_comm->MyPID() == 0 )
	{
		std::string this_prefix = "/"+M_prefix+"/";
		HDF5dense_exporter.WriteVectorDouble(this_prefix, "SingularValues", M_SingularValues);
	}

	if (usePrefix)
	{
		HDF5_exporter.Write("V"+M_prefix, *M_PODbasis);// Save solution basis V
		HDF5_exporter.Write("info", "num_basis_"+M_prefix, M_NumPODbases);
	}
	else
	{
		HDF5_exporter.Write("V", *M_PODbasis);// Save solution basis V
		HDF5_exporter.Write("info", "num_basis", M_NumPODbases);
	}

}
//=========================================================================
void ReductionVector::getPODbasis( boost::shared_ptr<Epetra_FECrsMatrix>& matrix )
{
	matrix = M_PODbasis;
}
//=========================================================================
} // end namespace LifeV
