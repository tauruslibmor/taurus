#include <taurus/Core/HyperReduction/ImportRhsSnapshots.hpp>
#include "Epetra_Vector.h"
#include <taurus/Core/Utilities/Timer.hpp>

namespace LifeV
{


ImportRhsSnapshots::ImportRhsSnapshots ( const mapPtr_Type& map_rows, const commPtr_Type& communicator, const int numParameters, const GetPot& dataFile, const std::string& dataFile_section ) :
	M_numParameters  ( numParameters ),
	M_datafile       ( dataFile ),
	M_map_rows		 ( map_rows ),
	M_comm           ( communicator ),
	M_datafileSection  ( dataFile_section )
{
}

void
ImportRhsSnapshots::readSnapshots( )
{

	if ( M_comm->MyPID() == 0 )
        std::cout << "\n Read " << M_numParameters << " " << M_datafileSection << " snapshots ... ";

	Timer myTimer(M_comm, 1.0, "s");
	myTimer.StartTimer();

	buildColumnMap();

	M_X.reset ( new Epetra_FECrsMatrix ( Copy, *M_map_rows->map(Unique), M_numParameters, false ) );

	std::string data_field = M_datafileSection+"/name_input_file";
	std::string fileName_importer =  M_datafile( data_field.c_str() , "SystemSnapshots")+".h5";

	EpetraExt::HDF5 HDF5(*M_comm);

	HDF5.Open( fileName_importer );

	data_field = M_datafileSection+"/prefix";
	std::string h5_rhs_prefix = M_datafile(data_field.c_str(), "F");
    
	Epetra_MultiVector* importedVector;

	for (int i = 0; i < M_numParameters; ++i )
	{

 		HDF5.Read(h5_rhs_prefix+"_"+EpetraExt::toString(i+1), *M_map_rows->map(Unique), importedVector, false );

		for ( int j = 0; j < M_map_rows->map(Unique)->NumMyElements(); ++j )
		{
			int index_row = M_map_rows->map(Unique)->GID(j);
			M_X->InsertGlobalValues( 1, &index_row, 1, &i, &(*(*importedVector)(0))[j] );
		}

	}

	M_X->GlobalAssemble(*M_map_column, *M_map_rows->map(Unique));

	HDF5.Close();

	double elapsedtime = myTimer.StopTimer();
	if ( M_comm->MyPID() == 0 )
			std::cout << "done in " << elapsedtime << " s";
}

void
ImportRhsSnapshots::buildColumnMap( )
{
	std::vector<int> column_indices;
	for ( int i = 0; i < M_numParameters; ++i )
		column_indices.push_back(i);

	int* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	//M_map_column.reset ( new MapEpetra ( -1, static_cast<int> (column_indices.size() ), pointerToDofs, M_comm ) );
	M_map_column.reset ( new Epetra_Map( (int)(-1), static_cast<int> (M_numParameters ), pointerToDofs, (int)(0), *M_comm ));

}


void
ImportRhsSnapshots::getMatrix ( FECrsPtr_Type& matrix )
{
	matrix.reset( new Epetra_FECrsMatrix ( *M_X ) );
	matrix->GlobalAssemble( *M_map_column, *M_map_rows->map(Unique) );
}

} // end namespace LifeV
