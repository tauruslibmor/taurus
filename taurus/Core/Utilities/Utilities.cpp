#include <taurus/Core/Utilities/Utilities.hpp>

namespace LifeV
{

void
RTRed_Utils::projectMatrix( const boost::shared_ptr<Epetra_FECrsMatrix>& A, const boost::shared_ptr<Epetra_FECrsMatrix>& V, const int num_bases,
                            const mapPtr_Type& map_rows, const mapPtr_Type& map_cols,  boost::shared_ptr<Epetra_FECrsMatrix>& out)
{
	boost::shared_ptr<Epetra_FECrsMatrix> tmp ( new Epetra_FECrsMatrix ( Copy, V->OperatorRangeMap(), num_bases ) );
	int errCode = EpetraExt::MatrixMatrix::Multiply ( *A, false, *V, false, *tmp, false );
	tmp->GlobalAssemble( V->OperatorDomainMap(), V->OperatorRangeMap() );
	out.reset ( new Epetra_FECrsMatrix ( Copy, V->OperatorDomainMap(), num_bases ) );
	errCode = EpetraExt::MatrixMatrix::Multiply ( *V, true, *tmp, false, *out, false );
	out->GlobalAssemble(V->OperatorDomainMap(),V->OperatorDomainMap());

//    matrixPtr_Type tmp ( new matrix_Type ( *map_rows, num_bases ) );
//    A->multiply(false, *V, false, *tmp, false);
//    tmp->globalAssemble( map_cols, map_rows );
//    out.reset ( new matrix_Type ( *map_cols, num_bases ) );
//    V->multiply(true, *tmp, false, *out, false);
//    out->globalAssemble();

}

void
RTRed_Utils::projectMatrix( const boost::shared_ptr<Epetra_FECrsMatrix>& A, const boost::shared_ptr<Epetra_FECrsMatrix>& V, const int num_bases,
                            boost::shared_ptr<Epetra_FECrsMatrix>& out)
{
	boost::shared_ptr<Epetra_FECrsMatrix> tmp ( new Epetra_FECrsMatrix ( Copy, V->OperatorRangeMap(), num_bases ) );
	int errCode = EpetraExt::MatrixMatrix::Multiply ( *A, false, *V, false, *tmp, false );
	tmp->GlobalAssemble( V->OperatorDomainMap(), V->OperatorRangeMap() );
	out.reset ( new Epetra_FECrsMatrix ( Copy, V->OperatorDomainMap(), num_bases ) );
	errCode = EpetraExt::MatrixMatrix::Multiply ( *V, true, *tmp, false, *out, false );
	out->GlobalAssemble(V->OperatorDomainMap(),V->OperatorDomainMap());

//    matrixPtr_Type tmp ( new matrix_Type ( *map_rows, num_bases ) );
//    A->multiply(false, *V, false, *tmp, false);
//    tmp->globalAssemble( map_cols, map_rows );
//    out.reset ( new matrix_Type ( *map_cols, num_bases ) );
//    V->multiply(true, *tmp, false, *out, false);
//    out->globalAssemble();

}

void
RTRed_Utils::transpose( boost::shared_ptr<Epetra_FECrsMatrix>& input_matrix, boost::shared_ptr<Epetra_CrsMatrix>& output_matrix)
{

	//boost::shared_ptr<Epetra_CrsMatrix> V_mat_transposed;
	Epetra_Map domain_map ( input_matrix->DomainMap() );

	EpetraExt::RowMatrix_Transpose transposer( &domain_map );
	output_matrix.reset( new Epetra_CrsMatrix ( *( dynamic_cast<Epetra_CrsMatrix*>(&(transposer(*input_matrix) ) ) ) ) );
	output_matrix->FillComplete( input_matrix->RangeMap(), input_matrix->DomainMap() );
}

void
RTRed_Utils::extractRow(const boost::shared_ptr<Epetra_CrsMatrix>& matrix, const int row_index, double* srcValues, int* colIndices)
{
	int numEntries;
	int numEntriesPerRow;
	int globalRow;

	numEntries = matrix->NumMyEntries(row_index);
	numEntriesPerRow = matrix->NumMyEntries(row_index);
	globalRow = matrix->OperatorRangeMap().GID64(row_index);

	matrix->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);

}

void
RTRed_Utils::FECrsMatrixToMultiVector(const boost::shared_ptr<Epetra_FECrsMatrix>& input_matrix, boost::shared_ptr<Epetra_MultiVector>& vector, int numcols)
{

	vector.reset( new Epetra_MultiVector (input_matrix->OperatorRangeMap(), numcols) );
	vector->PutScalar(0.0);

	int numEntries;
	int numEntriesPerRow;
	int globalRow;

	for (int l = 0; l < input_matrix->NumMyRows() ; l++ )
	{
		numEntries = input_matrix->NumMyEntries(l);
		numEntriesPerRow = input_matrix->NumMyEntries(l);
		globalRow = input_matrix->OperatorRangeMap().GID(l);
		double * srcValues = new double[numEntries];
		int * colIndices = new int[numEntries];

		input_matrix->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);

		for (int j = 0; j < input_matrix->NumMyEntries(l); j++ )
		{
			//vector->ReplaceGlobalValue(globalRow, colIndices[j], srcValues[j]);
			vector->ReplaceMyValue (l, colIndices[j], srcValues[j]);
		}

		delete [] srcValues;
		delete [] colIndices;
	}
}
    
} // namespace LifeV

