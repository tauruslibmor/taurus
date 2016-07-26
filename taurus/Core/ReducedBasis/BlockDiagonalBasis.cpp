#include <taurus/Core/ReducedBasis/BlockDiagonalBasis.hpp>
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
BlockDiagonalBasis::BlockDiagonalBasis (const commPtr_Type& communicator,
										const boost::shared_ptr<Epetra_FECrsMatrix>& V1,
										const int num_rows_V1,
										const boost::shared_ptr<Epetra_FECrsMatrix>& V2,
										const boost::shared_ptr<Epetra_Map>& target_RowMap ) :
	M_comm ( communicator ),
	M_map_rows ( target_RowMap ),
	M_V1 ( V1 ),
	M_V2 ( V2 ),
	M_num_cols_V ( 0 ),
	M_num_rows_V1 (num_rows_V1)
{}
//=========================================================================
void BlockDiagonalBasis::blockRow( )
{
	// mappa colonne
	int num_cols_V1 = M_V1->NumMyCols();
	int num_cols_V2 = M_V2->NumMyCols();
	M_num_cols_V = num_cols_V1+ num_cols_V2;

	// build column map for the basis V
	std::vector<int> column_indices;
	for ( int i = 0; i < M_num_cols_V; ++i )
		column_indices.push_back(i);

	int* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	boost::shared_ptr<Epetra_Map> map_column_unique ( new Epetra_Map( -1, static_cast<int> (column_indices.size() ), pointerToDofs, 0, *M_comm ) );

	M_V.reset ( new Epetra_FECrsMatrix ( Copy, *M_map_rows, M_num_cols_V, false ) );

	for (int l = 0; l < M_V1->NumMyRows() ; l++ )
	{
		int numEntries = M_V1->NumMyEntries(l);
		int numEntriesPerRow = M_V1->NumMyEntries(l);
		int globalRow = M_V1->OperatorRangeMap().GID(l);
		double * srcValues = new double[numEntries];
		int * colIndices = new int[numEntries];

		M_V1->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);

		for (int j = 0; j < numEntries; j++ )
		{
			M_V->InsertGlobalValues ( 1, &globalRow, 1, &colIndices[j], &srcValues[j] );
		}

		delete [] srcValues;
		delete [] colIndices;
	}

	for (int l = 0; l < M_V2->NumMyRows() ; l++ )
	{
		int numEntries = M_V2->NumMyEntries(l);
		int numEntriesPerRow = M_V2->NumMyEntries(l);
		int globalRow = M_V2->OperatorRangeMap().GID(l);
		double * srcValues = new double[numEntries];
		int * colIndices = new int[numEntries];

		M_V2->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);

		for (int j = 0; j < numEntries; j++ )
		{
			int row = globalRow;
			int col = colIndices[j]+num_cols_V1;
			M_V->InsertGlobalValues ( 1, &row, 1, &col, &srcValues[j] );
		}

		delete [] srcValues;
		delete [] colIndices;
	}

	M_V->GlobalAssemble(*map_column_unique, *M_map_rows);

	//std::string namem = "basisVBD.m";
	//EpetraExt::RowMatrixToMatlabFile ( namem.c_str(), *M_V );
}
//=========================================================================
void BlockDiagonalBasis::perform( )
{
	// mappa colonne
	int num_cols_V1 = M_V1->NumMyCols();
	int num_cols_V2 = M_V2->NumMyCols();
	M_num_cols_V = num_cols_V1+ num_cols_V2;

	// build column map for the basis V
	std::vector<int> column_indices;
	for ( int i = 0; i < M_num_cols_V; ++i )
		column_indices.push_back(i);

	int* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	boost::shared_ptr<Epetra_Map> map_column_unique ( new Epetra_Map( -1, static_cast<int> (column_indices.size() ), pointerToDofs, 0, *M_comm ) );

	M_V.reset ( new Epetra_FECrsMatrix ( Copy, *M_map_rows, M_num_cols_V, false ) );

	for (int l = 0; l < M_V1->NumMyRows() ; l++ )
	{
		int numEntries = M_V1->NumMyEntries(l);
		int numEntriesPerRow = M_V1->NumMyEntries(l);
		int globalRow = M_V1->OperatorRangeMap().GID(l);
		double * srcValues = new double[numEntries];
		int * colIndices = new int[numEntries];

		M_V1->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);

		for (int j = 0; j < numEntries; j++ )
		{
			M_V->InsertGlobalValues ( 1, &globalRow, 1, &colIndices[j], &srcValues[j] );
		}

		delete [] srcValues;
		delete [] colIndices;
	}

	for (int l = 0; l < M_V2->NumMyRows() ; l++ )
	{
		int numEntries = M_V2->NumMyEntries(l);
		int numEntriesPerRow = M_V2->NumMyEntries(l);
		int globalRow = M_V2->OperatorRangeMap().GID(l);
		double * srcValues = new double[numEntries];
		int * colIndices = new int[numEntries];

		M_V2->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);

		for (int j = 0; j < numEntries; j++ )
		{
			int row = globalRow+M_num_rows_V1;
			int col = colIndices[j]+num_cols_V1;
			M_V->InsertGlobalValues ( 1, &row, 1, &col, &srcValues[j] );
		}

		delete [] srcValues;
		delete [] colIndices;
	}

	M_V->GlobalAssemble(*map_column_unique, *M_map_rows);

	//std::string namem = "basisVBD.m";
	//EpetraExt::RowMatrixToMatlabFile ( namem.c_str(), *M_V );
}
//=========================================================================
void BlockDiagonalBasis::write( EpetraExt::HDF5& HDF5_exporter, const std::string& name )
{
	HDF5_exporter.Write(name, *M_V);// Save solution basis V
	HDF5_exporter.Write("info", "num_basis_"+name, M_num_cols_V);
}
//=========================================================================
void BlockDiagonalBasis::getBasis( boost::shared_ptr<Epetra_FECrsMatrix>& matrix )
{
	matrix = M_V;
}
//=========================================================================
} // end namespace LifeV
