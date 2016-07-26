#include <taurus/Core/HyperReduction/VectorizedToMatrices.hpp>

namespace LifeV
{


VectorizedToMatrices::VectorizedToMatrices ( ) :
		hasTranspose (false)
{
}

void
VectorizedToMatrices::transpose( )
{
	if(!hasTranspose)
	{
		Epetra_Map range_map ( M_vectorizedMatrix->RangeMap() );

		EpetraExt::RowMatrix_Transpose transposer( &range_map );
		M_vectorizedMatrix_transposed.reset( new Epetra_CrsMatrix ( *( dynamic_cast<Epetra_CrsMatrix*>(&(transposer(*M_vectorizedMatrix) ) ) ) ) );

		hasTranspose = true;
	}
}


void
VectorizedToMatrices::getFullMatrix ( const int& column_index, boost::shared_ptr<Epetra_FECrsMatrix>& matrix )
{
	matrix.reset( new Epetra_FECrsMatrix ( Copy, *M_map_matrix_rows, 150, false ) );

	// 1) traspongo matrice
	transpose();

	// 2) estraggo riga desiderata, indice colonna e' ora indice riga
	int numEntries = M_vectorizedMatrix_transposed->NumMyEntries(column_index);
	int numEntriesPerRow = M_vectorizedMatrix_transposed->NumMyEntries(column_index);
	long long int globalRow = M_vectorizedMatrix_transposed->OperatorRangeMap().GID64(column_index);
	double * srcValues = new double[numEntries];
	long long int * colIndices = new long long int[numEntries];

	M_vectorizedMatrix_transposed->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);

	// 3) inserisco nella matrice FE i valori estratti precedentemente
	int row_index, col_index;
	double value;

	for ( int i = 0; i < numEntries; ++i )
	{
		row_index = 0;
		col_index = 0;
		reshapeIndex(colIndices[i], row_index, col_index );
		value = srcValues[i];
//		std::cout << "ROW = " << row_index << ", COL = " << col_index << ", VALUE = " << value << "\n";
		matrix->InsertGlobalValues ( 1, &row_index, 1, &col_index, &value  );
	}

	matrix->GlobalAssemble(*M_map_matrix_rows,*M_map_matrix_rows);

	delete [] srcValues;
	delete [] colIndices;
}

void
VectorizedToMatrices::reshapeIndex(const long long int& index_in, int& index_row, int& index_col )
{
	int c = (int)(index_in/M_offset);
	int r = (int)(index_in) - (c*M_offset);

	index_row = r;
	index_col = c;
}

} // end namespace LifeV
