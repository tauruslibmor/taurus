#include <taurus/Core/HyperReduction/ImportMatrixSnapshots.hpp>


//#include "EpetraExt_Version.h"

#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"

#include "Epetra_Time.h"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "Trilinos_Util.h"
#include <taurus/Core/Utilities/Timer.hpp>



namespace LifeV
{


ImportMatrixSnapshots::ImportMatrixSnapshots ( const mapPtr_Type& map_rows, const commPtr_Type& communicator, const int numParameters, const GetPot& dataFile, const std::string& dataFile_section ) :
    M_numParameters  ( numParameters ),
    M_datafileSection  ( dataFile_section ),
	M_datafile       ( dataFile ),
	M_map_matrix_rows		 ( map_rows ),
    M_offset	     ( 0 ),
    M_comm           ( communicator )
{
}

void
ImportMatrixSnapshots::readSnapshots( )
{

    if ( M_comm->MyPID() == 0 )
        std::cout << "\n Read " << M_numParameters << " " << M_datafileSection << " snapshots ... ";
    
    Timer myTimer(M_comm, 1.0, "s");
    myTimer.StartTimer();
    
	buildColumnMap();

	std::string data_field = M_datafileSection+"/name_input_file";
	std::string fileName_importer =  M_datafile( data_field.c_str() , "SystemSnapshots")+".h5";

    EpetraExt::HDF5 HDF5(*M_comm);
	HDF5.Open( fileName_importer );

	data_field = M_datafileSection+"/prefix";
	std::string h5_matrix_prefix = M_datafile(data_field.c_str(), "A");

	Epetra_CrsMatrix* importedMatrix_map;

	HDF5.Read(h5_matrix_prefix+"_map", *M_map_matrix_rows->map(Unique), *M_map_matrix_rows->map(Unique), importedMatrix_map );

	transpose ( importedMatrix_map );

	buildRowMap ( M_transposed_matrix );

	M_X.reset ( new Epetra_FECrsMatrix ( Copy, *M_map_rows, M_numParameters, false ) );

	long long int sizeMatrix = (importedMatrix_map->OperatorRangeMap().MaxAllGID () + 1 );

	M_offset = (int)sizeMatrix;

	for (int i = 0; i < M_numParameters; ++i )
	{
		Epetra_CrsMatrix* importedMatrix;

		// Step1: leggere matrici in formato matriciale da HDF5
		HDF5.Read(h5_matrix_prefix+"_"+EpetraExt::toString(i+1), *M_map_matrix_rows->map(Unique), *M_map_matrix_rows->map(Unique), importedMatrix );


        transpose ( importedMatrix );
        
        // Step2: piazzare matrice vettorizzata nella i-sima colonna della matrice M_X
        
        int k = 0;
        int numEntries;
        int numEntriesPerRow;
        int globalRow;
        

		for (int l = 0; l< M_transposed_matrix->NumMyRows() ; l++ )
		{
			numEntries = M_transposed_matrix->NumMyEntries(l);
			numEntriesPerRow = M_transposed_matrix->NumMyEntries(l);
			globalRow = M_transposed_matrix->OperatorRangeMap().GID(l);
			double * srcValues = new double[numEntries];
			int * colIndices = new int[numEntries];

			M_transposed_matrix->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);

			for (int j = 0; j < M_transposed_matrix->NumMyEntries(l); j++ )
			{
				long long int rows = (long long int)(colIndices[j]) + sizeMatrix * ((long long int)(globalRow));//(*M_indexes_map_matrix_vectorized)[k];
				long long int cols = i;// colonna corrisponde a i-sima lettura
				double values = srcValues[j];

				M_X->InsertGlobalValues ( 1, &rows, 1, &cols, &values  );

				k = k + 1;
			}

			delete [] srcValues;
			delete [] colIndices;
		}
		
        M_transposed_matrix.reset();
        importedMatrix->Scale(0.0);
	}

	M_X->GlobalAssemble(*M_map_column, *M_map_rows);
	HDF5.Close();
    
    double elapsedtime = myTimer.StopTimer();
    if ( M_comm->MyPID() == 0 )
        std::cout << "done in " << elapsedtime << " s";
}

void
ImportMatrixSnapshots::buildColumnMap( )
{
	std::vector<long long int> column_indices;
	for ( long long int i = 0; i < ( long long int ) (M_numParameters); ++i )
		column_indices.push_back(i);

	long long int* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	//M_map_column.reset ( new Epetra_Map( (long long int)(-1), static_cast<int> (M_numParameters ), (long long int)(0), *M_comm ));
	M_map_column.reset ( new Epetra_Map( (long long int)(-1), static_cast<int> (M_numParameters ), pointerToDofs, (long long int)(0), *M_comm ));
}

void
ImportMatrixSnapshots::transpose( Epetra_CrsMatrix* input_matrix)
{
    //transposed_matrix.reset (new Epetra_FECrsMatrix (Copy, input_matrix->OperatorDomainMap(), input_matrix->OperatorRangeMap(), 0, false) );
	Epetra_Map range_map ( input_matrix->RangeMap() );
    EpetraExt::RowMatrix_Transpose transposer(&range_map);
    // *dynamic_cast<Epetra_CrsMatrix*> (& (*transposed_matrix) ) = dynamic_cast<Epetra_CrsMatrix&> (transposer (*input_matrix) );
    //transposed_matrix->FillComplete();

    M_transposed_matrix.reset( new Epetra_CrsMatrix ( *( dynamic_cast<Epetra_CrsMatrix*>(&(transposer(*input_matrix) ) ) ) ) );
    M_transposed_matrix->FillComplete();
}
    
void
ImportMatrixSnapshots::buildRowMap ( const boost::shared_ptr<Epetra_CrsMatrix>& matrix )
{
	// Mappa Vettorizzata: da implementare usando codice copy-paste qui sotto

    M_indexes_map_matrix_vectorized.reset ( new std::vector<long long int> ( ) );
	long long int sizeMatrix = (matrix->OperatorRangeMap().MaxAllGID () + 1 );

	int numEntriesPerRow;
	int numEntries;
	int globalRow;
	int k = 0;

	for (int i = 0; i< matrix->NumMyRows() ; i++ )
	{
		numEntries = matrix->NumMyEntries(i);
		numEntriesPerRow = matrix->NumMyEntries(i);
		globalRow = matrix->OperatorRangeMap().GID(i);
		double * srcValues = new double[numEntries];
		int * colIndices = new int[numEntries];

		matrix->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);

		for (int j = 0; j < matrix->NumMyEntries(i); j++ )
		{
			M_indexes_map_matrix_vectorized->push_back( (long long int)(colIndices[j]) + sizeMatrix * ((long long int)(globalRow)) );
		}

		delete [] srcValues;
		delete [] colIndices;

	}

	// Mappa matrice vettorizzata
	long long int* pointerToDofs (0);
	if ( M_indexes_map_matrix_vectorized->size() > 0)
		pointerToDofs = &((*M_indexes_map_matrix_vectorized)[0]);

    M_map_rows.reset ( new Epetra_Map ( (long long int)(-1), static_cast<int> (M_indexes_map_matrix_vectorized->size() ), pointerToDofs, (long long int)(0),
    		*M_comm ) );
}
    
void
ImportMatrixSnapshots::getMatrix ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix )
{
	matrix.reset( new Epetra_FECrsMatrix ( *M_X ) );

	//boost::shared_ptr<MapEpetra> map_rows ( new MapEpetra ( &M_X->OperatorRangeMap() ) );

	//boost::shared_ptr<MapEpetra> map_column ( new MapEpetra ( &M_X->OperatorDomainMap() ) );

	matrix->GlobalAssemble( M_X->OperatorDomainMap(), M_X->OperatorRangeMap() );
}
    
} // end namespace LifeV
