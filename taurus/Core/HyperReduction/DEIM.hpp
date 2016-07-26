//@HEADER
/*!
    @file
    @author D. Forti
    @date 10-11-2015
 */

#ifndef DEIM_H_
#define DEIM_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <taurus/Core/Utilities/DOF_Extractor.hpp>
#include <taurus/Core/Utilities/Timer.hpp>
#include "EpetraExt_Utils.h"
#include <EpetraExt_Transpose_RowMatrix.h>

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

#include "EpetraExt_HDF5.h"
#include "EpetraExt_Utils.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseSVD.h"

#include <iterator>
#include <chrono>

namespace LifeV
{

template <typename GIDtype>
class DEIM
{
public:

	typedef MapEpetra map_Type;
	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;

	typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<MatrixEpetra<Real>> matrixPtr_Type;

	DEIM ( const commPtr_Type& communicator, const boost::shared_ptr<Epetra_FECrsMatrix> matrix, const int numCols);

	DEIM ( const commPtr_Type& communicator, const Epetra_SerialDenseMatrix& A, const std::vector<GIDtype>& indices);

	~DEIM() {};

	void performDEIM();

	std::vector<GIDtype> getIndices() { return M_deim_GID; };

	void spy ( std::string const& fileName );

	Epetra_SerialDenseMatrix getInterpolationMatrix() { return M_interpolationMatrix; };

	void setInterpolationMatrix(const Epetra_SerialDenseMatrix& A) { M_interpolationMatrix = A; };

	void solveOnline(Epetra_SerialDenseVector& b, Epetra_SerialDenseVector& x);

	void extractEntriesVector(const vectorPtr_Type& vec, Epetra_SerialDenseVector& vec_DEIM);

	void extractEntriesMatrix(const matrixPtr_Type& matrix, const std::vector<int>& row_index, const std::vector<int>& col_index, Epetra_SerialDenseVector& vec_DEIM);


private:

	void transpose( boost::shared_ptr<Epetra_FECrsMatrix>& input_matrix, boost::shared_ptr<Epetra_CrsMatrix>& output_matrix);

	void transposeFE( boost::shared_ptr<Epetra_FECrsMatrix>& input_matrix, boost::shared_ptr<Epetra_FECrsMatrix>& output_matrix);

	void extractRow(const boost::shared_ptr<Epetra_CrsMatrix>& matrix, const int row_index, double* srcValues, GIDtype* colIndices);

	void maxIndexAbsVector(const double* srcValues, const int numEntries, const GIDtype* colIndices, double& max_value,
						   GIDtype& index_max_value, GIDtype& localColIndex, int& max_rank);

	void fillInterpolationMatrix(const boost::shared_ptr<Epetra_CrsMatrix>& V_mat_transposed,
								 const std::vector<int>& procs,  const std::vector<GIDtype>& deim_LID);

	void fillRhs(const boost::shared_ptr<Epetra_CrsMatrix>& V_mat_transposed, const std::vector<int>& procs,
			     const std::vector<GIDtype>& deim_LID, const int& column, const int& numEntries, Epetra_SerialDenseVector& b);

	void solveLinearSystem(Epetra_SerialDenseMatrix& A, Epetra_SerialDenseVector& b, Epetra_SerialDenseVector& x);

	void computeResidual(const boost::shared_ptr<Epetra_CrsMatrix>& V_mat_transposed, const Epetra_SerialDenseVector& x,
			   const int& numEntries, const int& column, double* residual);

	commPtr_Type M_comm;
	boost::shared_ptr<Epetra_FECrsMatrix> M_matrix;
	int M_numCols;
	std::vector<GIDtype> M_deim_GID;
	Epetra_SerialDenseMatrix M_interpolationMatrix;
};

// ==========================
// IMPLEMENTAZIONE
// ==========================

//=========================================================================
// Offline constructor
template <typename GIDtype>
DEIM<GIDtype>::DEIM ( const commPtr_Type& communicator, const boost::shared_ptr<Epetra_FECrsMatrix> matrix, const int numCols ) :
   M_comm  		( communicator ),
   M_matrix   		( matrix ),
   M_numCols       ( numCols )
{
}
//=========================================================================
// Online constructor
template <typename GIDtype>
DEIM<GIDtype>::DEIM ( const commPtr_Type& communicator, const Epetra_SerialDenseMatrix& A, const std::vector<GIDtype>& indices) :
   M_comm  		                ( communicator ),
   M_interpolationMatrix   		( A ),
   M_deim_GID                   ( indices )
{
}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::spy ( std::string const& fileName)
{
	std::string name = "", uti = " , ";

	Int  me = M_matrix->Comm().MyPID();
	std::ostringstream myStream;
	myStream << me;
	name = fileName + ".m";

	EpetraExt::RowMatrixToMatlabFile ( name.c_str(), *M_matrix );
}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::transpose( boost::shared_ptr<Epetra_FECrsMatrix>& input_matrix, boost::shared_ptr<Epetra_CrsMatrix>& output_matrix)
{
	Epetra_Map domain_map ( input_matrix->DomainMap() );

	EpetraExt::RowMatrix_Transpose transposer( &domain_map );
	output_matrix.reset( new Epetra_CrsMatrix ( *( dynamic_cast<Epetra_CrsMatrix*>(&(transposer(*input_matrix) ) ) ) ) );
	output_matrix->FillComplete( input_matrix->RangeMap(), domain_map);
}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::transposeFE( boost::shared_ptr<Epetra_FECrsMatrix>& input_matrix, boost::shared_ptr<Epetra_FECrsMatrix>& output_matrix)
{

	Timer myTimer(M_comm, 1.0, "s");

	std::vector<GIDtype> column_indices;
	for ( GIDtype i = 0; i < (GIDtype)(M_numCols); ++i )
		column_indices.push_back(i);

	GIDtype* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	Epetra_Map domain_map( static_cast<GIDtype> (column_indices.size() ), static_cast<GIDtype> (column_indices.size() ), pointerToDofs, (GIDtype)(0), *M_comm) ;

 	output_matrix.reset ( new Epetra_FECrsMatrix ( Copy, domain_map, input_matrix->RangeMap(), 0, false ) );

	int numEntries;
	int numEntriesPerRow;
	GIDtype globalRow;

	myTimer.StartTimer();
	for (GIDtype l = 0; l< (GIDtype)(input_matrix->NumMyRows()) ; l++ )
	{
 
		numEntries = input_matrix->NumMyEntries(l);
		numEntriesPerRow = input_matrix->NumMyEntries(l);
		globalRow = input_matrix->OperatorRangeMap().GID64(l);
		double * srcValues = new double[numEntries];
		GIDtype * colIndices = new GIDtype[numEntries];
 
		input_matrix->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);
 		for (int j = 0; j < input_matrix->NumMyEntries(l); j++ )
		{
			double values = srcValues[j];
			output_matrix->InsertGlobalValues ( colIndices[j], 1, &values, &globalRow );
		}
		//output_matrix->InsertGlobalValues ( input_matrix->NumMyEntries(l), colIndices, 1, &globalRow, &srcValues);

 		//InsertGlobalValues (int numRows, const long long *rows, int numCols, const long long *cols, const double *values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)

 		delete [] srcValues;
		delete [] colIndices;
	}
	myTimer.StopTimer("Loop");
 
	output_matrix->GlobalAssemble( input_matrix->RangeMap(), domain_map );
}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::extractRow(const boost::shared_ptr<Epetra_CrsMatrix>& matrix, const int row_index, double* srcValues, GIDtype* colIndices)
{
	int numEntries;
	int numEntriesPerRow;
	int globalRow;

	numEntries = matrix->NumMyEntries(row_index);
	numEntriesPerRow = matrix->NumMyEntries(row_index);
	globalRow = matrix->OperatorRangeMap().GID64(row_index);

	matrix->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);
}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::maxIndexAbsVector(const double* srcValues, const int numEntries, const GIDtype* colIndices, double& max_value, GIDtype& index_max_value, GIDtype& localColIndex, int& max_rank)
{
	// 1) ricerca locale
	max_value = -1.0;
	index_max_value = 0;
	localColIndex   = 0; // to handle empty bases

	for ( GIDtype i = 0; i < (GIDtype)(numEntries); ++i )
	{
		if ( std::abs(srcValues[i]) > max_value )
		{
			max_value = std::abs(srcValues[i]);
			index_max_value = colIndices[i];
			localColIndex = i;
		}
	}
	//std::cout << "\nProc " << M_comm->MyPID() << ": max_value = " << max_value << ", index_max_value = " << index_max_value << "\n";

	// 2) MPI_Reduce
	int myrank, root;

	struct {
		double val;
		int   rank;
		GIDtype   localCol;
	} in[1], out[1];
	root = 0;

	out[0].val = 0;
	out[0].rank = 0;
	out[0].localCol = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	in[0].val = max_value;
	in[0].rank = myrank;
	in[0].localCol = localColIndex;

	MPI_Reduce( in, out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, root, MPI_COMM_WORLD );

	MPI_Bcast(&out[0].val, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

	MPI_Bcast(&out[0].rank, 1, MPI_INT, root, MPI_COMM_WORLD);

	MPI_Bcast(&localColIndex, 1, MPI_LONG_LONG_INT, out[0].rank, MPI_COMM_WORLD);

	MPI_Bcast(&index_max_value, 1, MPI_LONG_LONG_INT, out[0].rank, MPI_COMM_WORLD);

	max_rank = out[0].rank;

}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::fillInterpolationMatrix(const boost::shared_ptr<Epetra_CrsMatrix>& V_mat_transposed, const std::vector<int>& procs, const std::vector<GIDtype>& deim_LID)
{

	// build serial dense matrix A
	M_interpolationMatrix.Reshape( M_deim_GID.size(), M_deim_GID.size() );
	M_interpolationMatrix.Scale(0.0);

	for ( int i_col = 0; i_col < M_deim_GID.size(); ++i_col )
	{
		int numEntries2 = V_mat_transposed->NumMyEntries(i_col);
		double* srcValues= new double[numEntries2];
		GIDtype* colIndices = new GIDtype[numEntries2];
		extractRow(V_mat_transposed, i_col, srcValues, colIndices);

		for ( int i_row = 0; i_row < M_deim_GID.size(); ++i_row )
		{
			double matrix_value = 0;
			if ( M_comm->MyPID() == procs[i_row] )
			{
				matrix_value = srcValues[deim_LID[i_row]];
			}

			MPI_Bcast(&matrix_value, 1, MPI_DOUBLE, procs[i_row], MPI_COMM_WORLD);

			M_interpolationMatrix(i_row,i_col) = matrix_value;
		}

		delete [] srcValues;
		delete [] colIndices;
	}

}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::fillRhs(const boost::shared_ptr<Epetra_CrsMatrix>& V_mat_transposed, const std::vector<int>& procs, const std::vector<GIDtype>& deim_LID,
		               const int& column, const int& numEntries, Epetra_SerialDenseVector& b)
{

	b.Size( M_deim_GID.size() );
	double* srcValues= new double[numEntries];
	GIDtype* colIndices = new GIDtype[numEntries];
	extractRow(V_mat_transposed, column, srcValues, colIndices);

	for ( int i_row = 0; i_row < M_deim_GID.size(); ++i_row )
	{
		double matrix_value = 0;
		if ( M_comm->MyPID() == procs[i_row] )
		{
			matrix_value = srcValues[deim_LID[i_row]];
		}

		MPI_Bcast(&matrix_value, 1, MPI_DOUBLE, procs[i_row], MPI_COMM_WORLD);

		b(i_row) = matrix_value;
	}

	delete [] srcValues;
	delete [] colIndices;

}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::solveLinearSystem(Epetra_SerialDenseMatrix& A, Epetra_SerialDenseVector& b, Epetra_SerialDenseVector& x)
{

	x.Size( M_deim_GID.size() );

	Epetra_SerialDenseSolver my_solver;

	my_solver.SetMatrix( A );
	my_solver.SetVectors( x , b );
	my_solver.Solve();

}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::solveOnline(Epetra_SerialDenseVector& b, Epetra_SerialDenseVector& x)
{
	x.Size( M_deim_GID.size() );
	Epetra_SerialDenseSolver my_solver;

	my_solver.SetMatrix( M_interpolationMatrix );
	my_solver.SetVectors( x , b );
    my_solver.Solve();

	/*
	// explicit calculation of the inverse
	x.Size( M_deim_GID.size() );
	Epetra_SerialDenseSVD my_solver;
	my_solver.SetMatrix( M_interpolationMatrix );
	my_solver.Invert();
	Epetra_SerialDenseMatrix* Inv;
	Inv = my_solver.InvertedMatrix();
	Inv->Multiply (false, b, x);
	*/
}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::computeResidual(const boost::shared_ptr<Epetra_CrsMatrix>& V_mat_transposed, const Epetra_SerialDenseVector& x,
							   const int& numEntries, const int& column, double* residual)
{

	// residual computation
	double * linear_combination = new double [numEntries];

	for ( int k = 0; k < numEntries; ++k )
	{
		residual[k] = 0;
		linear_combination[k] = 0;
	}

	for ( int i_col = 0; i_col < M_deim_GID.size(); ++i_col )
	{
		double* srcValues= new double[numEntries];
		GIDtype* colIndices = new GIDtype[numEntries];
		extractRow(V_mat_transposed, i_col, srcValues, colIndices);

		for ( int k = 0; k < numEntries; ++k )
		{
			linear_combination[k] += x(i_col) * srcValues[k];
		}

		delete [] srcValues;
		delete [] colIndices;
	}

	double* srcValues= new double[numEntries];
	GIDtype* colIndices = new GIDtype[numEntries];
	extractRow(V_mat_transposed, column, srcValues, colIndices);

	for ( int k = 0; k < numEntries; ++k )
	{
		residual[k] = srcValues[k] - linear_combination[k];
	}

	delete [] srcValues;
	delete [] colIndices;
	delete [] linear_combination;
}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::performDEIM()
{

	Timer myTimer(M_comm, 1.0, "s");

	if ( M_comm->MyPID() == 0 )
	{
		std::cout << "\n========================================";
		std::cout << "\nDEIM starts";
	}

	myTimer.StartGlobalTimer();

	std::vector<GIDtype> deim_LID;
	std::vector<int> procs;

	// Transpose matrix
	/*
	myTimer.StartTimer();
	boost::shared_ptr<Epetra_FECrsMatrix> V_mat_transposed;
	transposeFE( M_matrix, V_mat_transposed);
	myTimer.StopTimer("transposing the POD basis");
	*/

	myTimer.StartTimer();
	boost::shared_ptr<Epetra_CrsMatrix> V_mat_transposed;
	transpose( M_matrix, V_mat_transposed);
	myTimer.StopTimer("transposing the POD basis");

	int numEntries = V_mat_transposed->NumMyEntries(0);

	// First column index
	double* srcValues= new double[numEntries];
	GIDtype* colIndices = new GIDtype[numEntries];

	extractRow(V_mat_transposed, 0, srcValues, colIndices);

	double max_value;
	GIDtype index_max_value;
	GIDtype localColIndex;
	int max_rank;

	maxIndexAbsVector(srcValues, numEntries, colIndices, max_value, index_max_value, localColIndex, max_rank);
	M_deim_GID.push_back(index_max_value);
	deim_LID.push_back(localColIndex);
	procs.push_back(max_rank);

	delete [] srcValues;

	// delete [] colIndices;
	// Assumiamo che tutte le colonne abbiamo stesso colIndices, inutile cancellarlo, riusato piu' sotto

	if ( M_comm->MyPID() == 0 )
		std::cout << "GID Index " << 0 << " = " << M_deim_GID[0] << ", ";

	// Loop over the columns
	for ( int l = 1; l < M_numCols; ++l)
	{

		int numEntries = V_mat_transposed->NumMyEntries(l);

		// build serial dense matrix A
		fillInterpolationMatrix(V_mat_transposed, procs, deim_LID);

		// build serial dense rhs b
		Epetra_SerialDenseVector b;
		fillRhs(V_mat_transposed, procs, deim_LID, l, numEntries, b);

		// solve Ax = b
		Epetra_SerialDenseVector x;
		solveLinearSystem(M_interpolationMatrix, b, x);

		// residual computation
		double * residual = new double [numEntries];
		computeResidual(V_mat_transposed, x, numEntries, l, residual);

		// Find next index
		maxIndexAbsVector(residual, numEntries, colIndices, max_value, index_max_value, localColIndex, max_rank);

		M_deim_GID.push_back(index_max_value);
		deim_LID.push_back(localColIndex);
		procs.push_back(max_rank);

		delete [] residual;

		if ( M_comm->MyPID() == 0 )
			std::cout << "GID Index " << l << " = " << M_deim_GID[l] << ", ";
	}
	delete [] colIndices;

	// build interpolation matrix
	fillInterpolationMatrix(V_mat_transposed, procs, deim_LID);

	myTimer.StopGlobalTimer("running DEIM algorithm");
	if ( M_comm->MyPID() == 0 )
			std::cout << "========================================\n";

}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::extractEntriesVector(const vectorPtr_Type& vec, Epetra_SerialDenseVector& vec_DEIM)
{

	vec_DEIM.Size( M_deim_GID.size() );;

	int root = 0;

	for ( int k = 0; k < M_deim_GID.size(); ++k )
	{
		double entry_rhs = 0;
		int rank_test = -1;
		int rank_test_output;

		if ( vec->blockMap().LID( M_deim_GID[k] ) != -1 )
		{
			rank_test = M_comm->MyPID();

			entry_rhs = (*vec)[M_deim_GID[k]];
		}

		MPI_Reduce( &rank_test, &rank_test_output, 1, MPI_INT, MPI_MAX, root, MPI_COMM_WORLD );
		MPI_Bcast(&rank_test_output, 1, MPI_INT, root, MPI_COMM_WORLD);
		MPI_Bcast(&entry_rhs, 1, MPI_DOUBLE, rank_test_output, MPI_COMM_WORLD);
		vec_DEIM(k) = entry_rhs;
	}
}
//=========================================================================
template <typename GIDtype>
void
DEIM<GIDtype>::extractEntriesMatrix(const matrixPtr_Type& matrix, const std::vector<int>& row_index, const std::vector<int>& col_index, Epetra_SerialDenseVector& vec_DEIM)
{

	vec_DEIM.Size( row_index.size() );;

	int root = 0;

	for ( int k = 0; k < row_index.size(); ++k )
	{
		double entry_matrix = 0;
		int rank_test = -1;
		int rank_test_output;

		if ( matrix->matrixPtr()->RangeMap().LID( row_index[k] ) != -1 )
		{
			rank_test = M_comm->MyPID();
			int numEntries = matrix->matrixPtr()->NumGlobalEntries(row_index[k]);
			int numEntriesPerRow = matrix->matrixPtr()->NumGlobalEntries(row_index[k]);
			double * srcValues = new double[numEntries];
			int * colIndices = new int[numEntries];

			matrix->matrixPtr()->ExtractGlobalRowCopy( row_index[k], numEntriesPerRow, numEntries, srcValues, colIndices);

			for (int j = 0; j < numEntries; ++j )
			{
				if ( colIndices[j] == col_index[k] )
				{
					entry_matrix = srcValues[j];
				}
			}

			delete [] srcValues;
			delete [] colIndices;
		}

		MPI_Reduce( &rank_test, &rank_test_output, 1, MPI_INT, MPI_MAX, root, MPI_COMM_WORLD );
		MPI_Bcast(&rank_test_output, 1, MPI_INT, root, MPI_COMM_WORLD);
		MPI_Bcast(&entry_matrix, 1, MPI_DOUBLE, rank_test_output, MPI_COMM_WORLD);
		vec_DEIM(k) = entry_matrix;
	}
}
//=========================================================================
} // namespace LifeV

#endif //  DEIM_H_
