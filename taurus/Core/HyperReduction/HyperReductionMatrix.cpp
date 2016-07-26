#include <taurus/Core/HyperReduction/HyperReductionMatrix.hpp>
#include <taurus/Core/HyperReduction/DEIM.hpp>
#include <taurus/Core/HyperReduction/ImportMatrixSnapshots.hpp>
#include <taurus/Core/ReducedBasis/V_POD.hpp>
#include <taurus/Core/HyperReduction/VectorizedToMatrices.hpp>
#include <taurus/Core/Utilities/Utilities.hpp>


//#include "Epetra_Vector.h"

using namespace std;

namespace LifeV
{

//=========================================================================
// Constructor //
HyperReductionMatrix::HyperReductionMatrix (const GetPot& dataFile, const std::string& dataFile_section, const commPtr_Type& communicator, const mapPtr_Type& map_rows) :
	M_comm             ( communicator ),
	M_datafile         ( dataFile ),
	M_datafileSection  ( dataFile_section ),
	M_NumSnapshots	   ( 0 ),
	M_NumPODbases      ( 0 ),
	M_map_rows		   ( map_rows )
{
	std::string data_field = M_datafileSection+"/prefix";
	M_prefix = M_datafile(data_field.c_str(), "A");
}
//=========================================================================
void HyperReductionMatrix::perform( )
{

	EpetraExt::HDF5 HDF5_importer(*M_comm);

	HDF5_importer.Open( M_datafile ("importer/name_input_file_system", "SystemSnapshots")+".h5" );
	HDF5_importer.Read("info", "num_"+M_prefix+"Snapshots", M_NumSnapshots);
	HDF5_importer.Close();

	ImportMatrixSnapshots MATRIXreader (M_map_rows, M_comm, M_NumSnapshots, M_datafile, M_datafileSection );

	MATRIXreader.readSnapshots( );

	MATRIXreader.getMatrix( M_Snapshots );

	M_offset = MATRIXreader.getOffset();

	// POD
	if ( M_comm->MyPID() == 0 )
	{
		std::cout << "\n POD on " << M_datafileSection << " snapshots ...\n";
	}

	V_POD<long long int> podMAT ( M_comm, M_Snapshots, M_NumSnapshots );

	std::string data_field = M_datafileSection+"/nortol";
	podMAT.performPOD_vector ( M_PODbasis, M_datafile(data_field.c_str(), 1e-2) );
 
    M_SingularValues = podMAT.getSingularValues();

    M_NumPODbases  = podMAT.getNumCols();

    M_Snapshots.reset();// clear memory

	// DEIM
	DEIM<long long int> DEIM_matrix (M_comm, M_PODbasis, M_NumPODbases );

    DEIM_matrix.performDEIM( );

	M_DEIM_Phi = DEIM_matrix.getInterpolationMatrix();

	M_deim_GID = DEIM_matrix.getIndices();

	for ( int k = 0; k < M_deim_GID.size(); ++k )
	{
		int c = (int)(M_deim_GID[k]/M_offset);
		int r = M_deim_GID[k] - (c)*M_offset;

		M_deim_row_indices.push_back(r);
		M_deim_col_indices.push_back(c);

		M_deim_dofs_GID.push_back(c);
		M_deim_dofs_GID.push_back(r);
	}

}
//=========================================================================
void HyperReductionMatrix::write( DenseHDF5& HDF5dense_exporter )
{
	if ( M_comm->MyPID() == 0 )
	{
		// Hyper-reduction HDF5 file manager
		//DenseHDF5 HDF5dense_exporter;
		//HDF5dense_exporter.Open( M_datafile("rom_exporter/hyper_output_file", "Hyper_ROM")+".h5" );

		std::string this_prefix = "/"+M_prefix+"/";
		// Write
		HDF5dense_exporter.WriteVectorDouble(this_prefix, "SingularValues", M_SingularValues);
		HDF5dense_exporter.WriteEpetraSerialDenseMatrix(this_prefix, "PHI", M_DEIM_Phi);
		HDF5dense_exporter.WriteVectorLongInt(this_prefix, "indices_vec", M_deim_GID);
		HDF5dense_exporter.WriteVectorInt(this_prefix, "row_indices", M_deim_row_indices);
		HDF5dense_exporter.WriteVectorInt(this_prefix, "col_indices", M_deim_col_indices);
		HDF5dense_exporter.WriteVectorInt(this_prefix, "dofs", M_deim_dofs_GID);
		HDF5dense_exporter.WriteIntValue(this_prefix, "size_matrix", M_offset);

		//HDF5dense_exporter.Close();
	}

}
//=========================================================================
void HyperReductionMatrix::getPODbasis ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix )
{
	matrix.reset( new Epetra_FECrsMatrix ( *M_PODbasis ) );
	matrix->GlobalAssemble( M_PODbasis->OperatorDomainMap(), M_PODbasis->OperatorRangeMap() );
}
//=========================================================================
void HyperReductionMatrix::project(const boost::shared_ptr<Epetra_FECrsMatrix>& V, const int numColsV, EpetraExt::HDF5& HDF5_exporter)
{

    Timer myTimer(M_comm, 1.0, "s");
    
	VectorizedToMatrices reshaper;
	reshaper.setVectorizedMatrix(M_PODbasis);
	reshaper.setMatrixRowMap( M_map_rows->map(Unique) );
	reshaper.setOffset(M_offset);

 	boost::shared_ptr<Epetra_FECrsMatrix> AN;
	boost::shared_ptr<Epetra_FECrsMatrix> FE_matrix;

	RTRed_Utils utils;

	for ( int i = 1; i <= M_NumPODbases; ++i )
	{
        myTimer.StartTimer();
		reshaper.getFullMatrix(i-1, FE_matrix);
		utils.projectMatrix( FE_matrix, V, numColsV, AN );
        double elaps1 = myTimer.StopTimer();
        //std::cout << *FE_matrix;
        //AN->GlobalAssemble();

        // sum contributes distributed over all processors and send it to Proc 0
        // see http://trilinos.sandia.gov/Trilinos10.6Tutorial.pdf page 22
        int NumMyElements_target;
        if( M_comm->MyPID() == 0 )
        	NumMyElements_target = numColsV;
        else
        	NumMyElements_target = 0;

        Epetra_Map TargetMap(-1,NumMyElements_target,0,*M_comm);
        Epetra_Export Exporter(AN->OperatorDomainMap(),TargetMap);
        boost::shared_ptr<Epetra_FECrsMatrix> AN_target ( new Epetra_FECrsMatrix ( Copy, TargetMap, numColsV ) );
        AN_target->Export(*AN,Exporter,Add);
        AN_target->GlobalAssemble(TargetMap,TargetMap);

        myTimer.StartTimer();
        HDF5_exporter.Write(M_prefix+"N_"+EpetraExt::toString(i), *AN_target);
        double elaps2 = myTimer.StopTimer();

        if (M_comm->MyPID()==0)
        {
            std::streamsize ss = cout.precision();
            cout << std::setprecision(2) << std::scientific;
            cout << "Matrix " << M_prefix+"N_" << i << " computed in " << elaps1
            << " s and exported in " <<  elaps2 << " s\n";
            cout.unsetf(ios_base::fixed);
            cout.precision(ss);
        }
	}
    
	if ( M_comm->MyPID() == 0 )
	{
		std::cout << "\n============================\n";
	}

	// close HDF5 file manager
	//HDF5_exporter.Close();
}
//=========================================================================
void HyperReductionMatrix::spySnapshots( )
{

}
//=========================================================================
void HyperReductionMatrix::spyPOD( )
{

}
//=========================================================================
} // end namespace LifeV
