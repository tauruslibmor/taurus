//@HEADER
/*!
    @file
    @author D. Forti
    @date 10-11-2015
*/

#ifndef V_POD_H_
#define V_POD_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/filter/GetPot.hpp>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Epetra_LAPACK.h"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

#include <taurus/Core/Utilities/Timer.hpp>

namespace LifeV
{

template <typename GIDtype>
class V_POD
{
public:

	typedef RegionMesh< LinearTetra > mesh_Type;
	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;
	typedef boost::shared_ptr<Epetra_Map> mapPtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<BCHandler> bcPtr_Type;
	typedef boost::shared_ptr<MatrixEpetra<Real>> matrixPtr_Type;

	V_POD ( const commPtr_Type& communicator, const matrixPtr_Type& matrix, const int numParameters);

	V_POD ( const commPtr_Type& communicator, const boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const int numParameters);

    ~V_POD() {};

    void setWeightMatrix ( boost::shared_ptr<Epetra_FECrsMatrix> matrix );

    void performPOD_vector ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const Real nortol, const Real& LOBPCG_tol, const int& LOBPCG_maxRestarts );

    void performPOD_vector ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const Real nortol );

    void performPOD_sol ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const Real nortol, const Real& LOBPCG_tol, const int& LOBPCG_maxRestarts );
    
    void performPOD_sol ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const Real nortol);

    void performPOD_matrix ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const Real nortol, const Real& LOBPCG_tol, const int& LOBPCG_maxRestarts );

    int getNumCols() const { return M_indexNortol; };

    std::vector<double> getSingularValues() const { return M_singularValues; };

private:

    void checkTolerance ( const Real nortol );

    void initializeV ( const int num_cols );

    void FECrs_to_SerialDense(const boost::shared_ptr<Epetra_FECrsMatrix>& A, Teuchos::SerialDenseMatrix<int, double>& A_dense );

    void GESVD_wrapper(const Teuchos::SerialDenseMatrix<int,double>& A_dense, Teuchos::SerialDenseMatrix<int, double>& eigenvectors, std::vector<double>& eigenvalues );
    
    void setIndexNortol(const std::vector<double>& eigenvalues );


    commPtr_Type M_comm;
    Real M_nortol;
    mapPtr_Type M_map_column_repeated;
    mapPtr_Type M_map_column_unique;
    boost::shared_ptr<Epetra_FECrsMatrix> M_matrix;
    boost::shared_ptr<Epetra_FECrsMatrix> M_weightMatrix;
    bool M_useWeightMatrix;
    int M_numParameters;
    int M_numBlocks;
    int M_nev;
    int M_indexNortol;
    Real M_LOBPCG_tol;
    int M_LOBPCG_maxRestarts;
    std::vector<double> M_singularValues;// updated only for Lapack version
};

// ==========================
// IMPLEMENTAZIONE
// ==========================

//=========================================================================================
template <typename GIDtype>
V_POD<GIDtype>::V_POD ( const commPtr_Type& communicator, const matrixPtr_Type& matrix, const int numParameters ) :
	M_comm  		( communicator ),
	M_matrix   		( matrix->matrixPtr() ),
	M_numParameters ( numParameters ),
	M_indexNortol   ( 0 ),
	M_useWeightMatrix (false)
{
	M_singularValues.resize( M_numParameters );
}
//=========================================================================================
template <typename GIDtype>
V_POD<GIDtype>::V_POD ( const commPtr_Type& communicator, const boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const int numParameters ) :
	M_comm  		( communicator ),
	M_matrix   		( matrix ),
	M_numParameters ( numParameters ),
	M_indexNortol   ( 0 ),
	M_useWeightMatrix (false)
{
	M_singularValues.resize( M_numParameters );
}
//=========================================================================================
template <typename GIDtype>
void
V_POD<GIDtype>::checkTolerance (  const Real nortol )
{
	M_nortol = nortol;

	if ( M_nortol < 1.0 )
	{
		M_nev = M_numParameters;
		M_numBlocks = M_numParameters;
	}
	else
	{
		M_nev = (int)M_nortol;
		M_numBlocks = M_nev + 1;
	}
}
//=========================================================================================
template <typename GIDtype>
void
V_POD<GIDtype>::setWeightMatrix ( boost::shared_ptr<Epetra_FECrsMatrix> matrix )
{
	M_weightMatrix = matrix;
	M_useWeightMatrix = true;
}
//=========================================================================================
template <typename GIDtype>
void
V_POD<GIDtype>::setIndexNortol(  const std::vector<double>& eigenvalues )
{
    Real denominator = 0.0;
    Real numerator = 0.0;
    int  numev = eigenvalues.size();
    
    for (int i=0; i<numev; i++)
    {
        denominator += eigenvalues[i];
    }
    
    if (M_nortol >= 1.0)
    {
        M_indexNortol = M_nortol;
    }
    else
    {
        M_indexNortol = -1;
        int k = 0;
        Real check = 0.0;
        Real test = 1.0 - ( M_nortol*M_nortol );
        
        while ( M_indexNortol < 0 && k <= numev )
        {
            numerator = 0.0;
            
            for ( int i = 0; i < k; ++i)
            {
                numerator += eigenvalues[i];
            }
            
            check = numerator/denominator;
            
            if ( check >= test )
            {
                M_indexNortol = k;
            }
            else
            {
                ++k;
            }
        }
    }

    if (M_indexNortol<0)
    	M_indexNortol = 0;


    if ( M_comm->MyPID() == 0 )
        std::cout << "\n Number of POD modes retained: " << M_indexNortol << "\n";
}
//=========================================================================================
template <typename GIDtype>
void
V_POD<GIDtype>::FECrs_to_SerialDense(const boost::shared_ptr<Epetra_FECrsMatrix>& A, Teuchos::SerialDenseMatrix<int, double>& A_dense )
{
	for (int i = 0; i< A->NumMyRows() ; i++ )
	{
		int numEntries = A->NumMyEntries(i);
		int numEntriesPerRow = A->NumMyEntries(i);
		int globalRow = A->OperatorRangeMap().GID64(i); //ATTENZIONE!!!

		double * srcValues = new double[numEntries];
		GIDtype * colIndices = new GIDtype[numEntries];

		A->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);

		for (int j = 0; j < A->NumMyEntries(i); j++ )
		{
			A_dense((int)i, (int)colIndices[j]) = srcValues[j];
		}

		delete [] srcValues;
		delete [] colIndices;
	}
}
//=========================================================================================
template <typename GIDtype>
void
V_POD<GIDtype>::GESVD_wrapper(const Teuchos::SerialDenseMatrix<int,double>& A_dense, Teuchos::SerialDenseMatrix<int, double>& eigenvectors, std::vector<double>& eigenvalues )
{
	// for a description of GESVD inputs see
	// http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing.html#ga84fdf22a62b12ff364621e4713ce02f2
	// and
	// https://trilinos.org/docs/dev/packages/epetra/doc/html/classEpetra__LAPACK.html
	// See also https://trilinos.org/docs/dev/packages/teuchos/doc/html/LAPACK_2cxx_main_8cpp-example.html#_a0

	Epetra_LAPACK lapack;
	int info, lwork;
	double Vt[ 1 ];
	lwork = EPETRA_MAX( 3*M_numParameters + M_numParameters, 5*M_numParameters );
	std::vector<double> work( lwork );

	Timer myTimer(M_comm, 1.0, "s");
	myTimer.setSerial();
	myTimer.StartTimer();
	char JOBU  = 'A';
	char JOBVT = 'N';
	lapack.GESVD( JOBU, JOBVT, M_numParameters, M_numParameters, A_dense.values(), A_dense.stride(), &eigenvalues[0], eigenvectors.values(), M_numParameters, Vt, 1, &work[0], &lwork, &info );


	std::cout << "\n==================================" << std::endl;
	myTimer.StopTimer("solving Eigen-Problem with LAPACK GESVD");

	std::cout << "eigenvalues: \n" << std::endl;
	for (int i = 0; i < eigenvalues.size() ; i++ )
	{
		std::cout << "   " << eigenvalues[i]  << std::endl;
	}

}
//=========================================================================================
//////////////////////////////// LOBPCG VERSION ////////////////////////////////
template <typename GIDtype>
void
V_POD<GIDtype>::performPOD_sol ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const Real nortol, const Real& LOBPCG_tol, const int& LOBPCG_maxRestarts  )
{
	checkTolerance(nortol);

	int blockSize = M_numParameters;
	int maxRestarts = LOBPCG_maxRestarts;
	int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
	Teuchos::LAPACK<int,double> lapack;
	double tol = LOBPCG_tol;
	std::string which = "LM";

	Teuchos::ParameterList MyPL;
	MyPL.set( "Verbosity", verbosity );
	MyPL.set( "Which", which );
	MyPL.set( "Block Size", blockSize );
	MyPL.set( "Num Blocks", M_numBlocks );
	MyPL.set( "Maximum Restarts", maxRestarts );
	MyPL.set( "Convergence Tolerance", tol );
	typedef Anasazi::MultiVec<double> MV;
	typedef Anasazi::Operator<double> OP;


	// Create an Anasazi::EpetraMultiVec for an initial vector to start the solver.
	// Note: This needs to have the same number of columns as the blocksize.
	Teuchos::RCP<Anasazi::EpetraMultiVec> ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec(M_matrix->OperatorDomainMap(), blockSize) );
	ivec->MvRandom();


	// Call the constructor for the (A^T*A) operator
	Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix( *M_matrix ) );
	Teuchos::RCP<Anasazi::EpetraSymOp> Amat = Teuchos::rcp( new Anasazi::EpetraSymOp( A ) );


	Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
			Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(Amat, ivec) );


	// Inform the eigenproblem that the matrix A is symmetric
	MyProblem->setHermitian(true);

	// Set the number of eigenvalues requested and the blocksize the solver should use
	MyProblem->setNEV( M_nev );

	// Inform the eigenproblem that you are finished passing it information
	bool boolret = MyProblem->setProblem();
	if (!boolret) {
		throw "Anasazi::BasicEigenproblem::setProblem() returned with error.";
	}

	// Initialize the Block Arnoldi solver
	Anasazi::LOBPCGSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

	// Solve the problem to the specified tolerances or length
	Anasazi::ReturnType returnCode = MySolverMgr.solve();
	if (returnCode != Anasazi::Converged && M_comm->MyPID() ==0) {
		std::cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << std::endl;
	}

	// Get the eigenvalues and eigenvectors from the eigenproblem
	Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
	std::vector<Anasazi::Value<double> > evals = sol.Evals;
	int numev = sol.numVecs;

	const double one = 1.0;
	const double zero = 0.0;

	//
	// Compute singular values which are the square root of the eigenvalues
	//

	if (M_comm->MyPID()==0)
	{
		std::cout<<"------------------------------------------------------"<<std::endl;
		std::cout<<"Computed Singular Values: "<<std::endl;
		std::cout<<"------------------------------------------------------"<<std::endl;
	}

	//
	// Check on nortol
	//

	Real denominator = 0.0;
	Real numerator = 0.0;

	for (int i=0; i<numev; i++)
	{
		denominator += evals[i].realpart;
	}


	if (M_nortol >= 1.0)
	{ M_indexNortol = M_nortol;
	}
	else
	{
		M_indexNortol = -1;
		int k = 0;
		Real check = 0.0;
		Real test = 1.0 - ( M_nortol*M_nortol );

		while ( M_indexNortol < 0 && k <= numev )
		{
			numerator = 0.0;

			for ( int i = 0; i < k; ++i)
			{
				numerator += evals[i].realpart;
			}

			check = numerator/denominator;

			if ( check >= test )
			{
				M_indexNortol = k;
			}
			else
			{
				++k;
			}
		}
	}

	for (int i=0; i<M_indexNortol; i++)
	{
		evals[i].realpart = Teuchos::ScalarTraits<double>::squareroot( evals[i].realpart );
	}
	//
	// Compute left singular vectors : u = Av/sigma
	//
	std::vector<double> tempnrm(M_indexNortol), directnrm(M_indexNortol);
	std::vector<int> index(M_indexNortol);
	for (int i=0; i<M_indexNortol; i++)
	{
		index[i] = i;
	}
	Anasazi::EpetraMultiVec Av(M_matrix->OperatorRangeMap(), M_indexNortol), u(M_matrix->OperatorRangeMap(), M_indexNortol);
	Anasazi::EpetraMultiVec* evecs = dynamic_cast<Anasazi::EpetraMultiVec* >(sol.Evecs->CloneViewNonConst( index ));

	Teuchos::SerialDenseMatrix<int,double> S(M_indexNortol,M_indexNortol);
	A->Apply( *evecs, Av );
	Av.MvNorm( tempnrm );
	for (int i=0; i<M_indexNortol; i++)
	{
		S(i,i) = one/tempnrm[i];
	};
	u.MvTimesMatAddMv( one, Av, S, zero );

	//
	// Compute direct residuals : || Av - sigma*u ||
	//

	for (int i=0; i<M_indexNortol; i++)
	{
		S(i,i) = evals[i].realpart;
	}

	Av.MvTimesMatAddMv( -one, u, S, one );
	Av.MvNorm( directnrm );
	if (M_comm->MyPID()==0)
	{
		std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
		std::cout<<std::setw(16)<<"Singular Value"
				<<std::setw(20)<<"Direct Residual"
				<<std::endl;
		std::cout<<"------------------------------------------------------"<<std::endl;
		for (int i=0; i<numev; i++)
		{
			std::cout<<std::setw(16)<<evals[i].realpart
					<<std::setw(20)<<directnrm[i]
					                           <<std::endl;
		}
		std::cout<<"------------------------------------------------------"<<std::endl;
	}

	initializeV(M_indexNortol);

	/*
	std::vector<GIDtype> column_indices;
	for ( GIDtype i = 0; i < (GIDtype)(M_indexNortol); ++i )
		column_indices.push_back(i);

	GIDtype* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	// ATTENZIONE!!!!!
	//	M_map_column_unique.reset ( new Epetra_Map( (GIDtype)(-1), static_cast<int> (column_indices.size() ), (GIDtype)(0), *M_comm) ); // prima
	M_map_column_unique.reset ( new Epetra_Map( (GIDtype)(-1), static_cast<int> (column_indices.size() ), pointerToDofs, (GIDtype)(0), *M_comm) ); // adesso

    */

	matrix.reset ( new Epetra_FECrsMatrix (Copy, M_matrix->OperatorRangeMap(), M_indexNortol ) );

	int * FirstPointInElementList1 = NULL;

	FirstPointInElementList1 = u.Map().FirstPointInElementList();

	double **A_Pointers = u.Pointers();

	for ( GIDtype i = 0; i < (GIDtype)(u.Map().NumMyElements()); i++ )
	{
		for ( GIDtype ii = 0; ii < (GIDtype)(u.Map().ElementSize(i)); ii++ )
		{
			GIDtype iii;
			iii = (GIDtype)(FirstPointInElementList1[i])+ii;
			for ( GIDtype j = 0; j < M_indexNortol; j++ )
			{
				GIDtype irow = (GIDtype)(u.Map().GID(i));
				GIDtype icol = (GIDtype)(j);

				matrix->InsertGlobalValues ( 1, &irow, 1, &icol, &A_Pointers[j][iii] );
			}
		}
	}

	matrix->GlobalAssemble(*M_map_column_unique, M_matrix->OperatorRangeMap());

}
//=========================================================================================
//////////////////////////////// LAPACK VERSION ////////////////////////////////
template <typename GIDtype>
void
V_POD<GIDtype>::performPOD_sol ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const Real nortol)
{
	Timer myTimer(M_comm, 1.0, "s");
	myTimer.StartGlobalTimer();

    checkTolerance(nortol);
    
    // A = M_matrix^T * M_matrix; in futuro A = M_matrix^T * matrice_norma *  M_matrix;
    boost::shared_ptr<Epetra_FECrsMatrix> A (new Epetra_FECrsMatrix(Copy, M_matrix->OperatorDomainMap(), 100));
    A->Scale(0.0);
    EpetraExt::MatrixMatrix::Multiply (*M_matrix, true, *M_matrix, false, *A, false);
    A->GlobalAssemble(M_matrix->OperatorDomainMap(),M_matrix->OperatorDomainMap());
    
    Teuchos::SerialDenseMatrix<int, double> eigenvectors( M_numParameters , M_numParameters );
    std::vector<double> eigenvalues(M_numParameters);
    
    // Processore 0 di sicuro possiede tutta la matrice distribuita, solo lui fa EIGENSOLVE
    int root = 0;
    
    if ( M_comm->MyPID() == root )
    {
        
        // convert A into a SerialDenseMatrix
        Teuchos::SerialDenseMatrix<int, double> A_dense( M_numParameters , M_numParameters );
        FECrs_to_SerialDense(A, A_dense );
        
        // run SVD on A_dense
        GESVD_wrapper(A_dense, eigenvectors, eigenvalues );
    }
    
    M_comm->Broadcast( eigenvectors.values(), eigenvectors.numRows()*eigenvectors.numCols(), root );
    M_comm->Broadcast( &eigenvalues[0], M_numParameters, root );
    
    // compute singular values
    //std::vector<double> M_singularValues(M_numParameters);
    for (int i = 0; i < eigenvalues.size() ; i++ )
    {
    	M_singularValues[i] = std::sqrt( eigenvalues[i] );
    }
    
    if ( M_comm->MyPID() == 0 )
    {
        std::cout << "\n==================================\n" << std::endl;
        std::cout << "singular values: \n" << std::endl;
        for (int i = 0; i < eigenvalues.size() ; i++ )
        {
            std::cout << "   " << M_singularValues[i]  << std::endl;
        }
    }
    
    // Check on nortol
    setIndexNortol( eigenvalues );
    
    initializeV(M_indexNortol);
    
    boost::shared_ptr<Epetra_FECrsMatrix> tmp (new Epetra_FECrsMatrix(Copy, M_matrix->OperatorDomainMap(), M_indexNortol));
    

    for ( GIDtype i = 0; i < (GIDtype) (M_matrix->OperatorDomainMap().NumMyElements()); i++ )
    {
    	for ( GIDtype j = 0; j < (GIDtype)(M_map_column_unique->NumMyElements()); j++ )
    	{
    		double value = eigenvectors( (int)(i),(int)(j) ) / M_singularValues[j]; //decommentare
    		tmp->InsertGlobalValues ( 1, &i, 1, &j, &value );
    	}
    }

    tmp->GlobalAssemble(*M_map_column_unique, M_matrix->OperatorDomainMap());
    
    // fill basis V
    matrix.reset ( new Epetra_FECrsMatrix (Copy, M_matrix->OperatorRangeMap(), M_indexNortol ) );
    
    EpetraExt::MatrixMatrix::Multiply ( *M_matrix, false, *tmp, false, *matrix, false );
    
    matrix->GlobalAssemble(*M_map_column_unique, M_matrix->OperatorRangeMap());

    myTimer.StopGlobalTimer("performing POD");
        
}
//=========================================================================================
//////////////////////////////// LOBPCG VERSION ////////////////////////////////
template <typename GIDtype>
void
V_POD<GIDtype>::performPOD_vector ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const Real nortol, const Real& LOBPCG_tol, const int& LOBPCG_maxRestarts  )
{
	checkTolerance(nortol);

	int blockSize = M_numParameters;
	int maxRestarts = LOBPCG_maxRestarts;
	int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
	Teuchos::LAPACK<int,double> lapack;
	double tol = LOBPCG_tol;
	std::string which = "LM";

	Teuchos::ParameterList MyPL;
	MyPL.set( "Verbosity", verbosity );
	MyPL.set( "Which", which );
	MyPL.set( "Block Size", blockSize );
	MyPL.set( "Num Blocks", M_numBlocks );
	MyPL.set( "Maximum Restarts", maxRestarts );
	MyPL.set( "Convergence Tolerance", tol );
	typedef Anasazi::MultiVec<double> MV;
	typedef Anasazi::Operator<double> OP;


	// Create an Anasazi::EpetraMultiVec for an initial vector to start the solver.
	// Note: This needs to have the same number of columns as the blocksize.
	Teuchos::RCP<Anasazi::EpetraMultiVec> ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec(M_matrix->OperatorDomainMap(), blockSize) );
	ivec->MvRandom();


	// Call the constructor for the (A^T*A) operator
	Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix( *M_matrix ) );
	Teuchos::RCP<Anasazi::EpetraSymOp> Amat = Teuchos::rcp( new Anasazi::EpetraSymOp( A ) );


	Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
			Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(Amat, ivec) );


	// Inform the eigenproblem that the matrix A is symmetric
	MyProblem->setHermitian(true);

	// Set the number of eigenvalues requested and the blocksize the solver should use
	MyProblem->setNEV( M_nev );

	// Inform the eigenproblem that you are finished passing it information
	bool boolret = MyProblem->setProblem();
	if (!boolret) {
		throw "Anasazi::BasicEigenproblem::setProblem() returned with error.";
	}

	// Initialize the Block Arnoldi solver
	Anasazi::LOBPCGSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

	// Solve the problem to the specified tolerances or length
	Anasazi::ReturnType returnCode = MySolverMgr.solve();
	if (returnCode != Anasazi::Converged && M_comm->MyPID() ==0) {
		std::cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << std::endl;
	}

	// Get the eigenvalues and eigenvectors from the eigenproblem
	Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
	std::vector<Anasazi::Value<double> > evals = sol.Evals;
	int numev = sol.numVecs;

	const double one = 1.0;
	const double zero = 0.0;

	//
	// Compute singular values which are the square root of the eigenvalues
	//

	if (M_comm->MyPID()==0)
	{
		std::cout<<"------------------------------------------------------"<<std::endl;
		std::cout<<"Computed Singular Values: "<<std::endl;
		std::cout<<"------------------------------------------------------"<<std::endl;
	}

	//
	// Check on nortol
	//

	Real denominator = 0.0;
	Real numerator = 0.0;

	for (int i=0; i<numev; i++)
	{
		denominator += evals[i].realpart;
	}


	if (M_nortol >= 1.0)
	{ M_indexNortol = M_nortol;
	}
	else
	{
		M_indexNortol = -1;
		int k = 0;
		Real check = 0.0;
		Real test = 1.0 - ( M_nortol*M_nortol );

		while ( M_indexNortol < 0 && k <= numev )
		{
			numerator = 0.0;

			for ( int i = 0; i < k; ++i)
			{
				numerator += evals[i].realpart;
			}

			check = numerator/denominator;

			if ( check >= test )
			{
				M_indexNortol = k;
			}
			else
			{
				++k;
			}
		}
	}

	for (int i=0; i<M_indexNortol; i++)
	{
		evals[i].realpart = Teuchos::ScalarTraits<double>::squareroot( evals[i].realpart );
	}
	//
	// Compute left singular vectors : u = Av/sigma
	//
	std::vector<double> tempnrm(M_indexNortol), directnrm(M_indexNortol);
	std::vector<int> index(M_indexNortol);
	for (int i=0; i<M_indexNortol; i++)
	{
		index[i] = i;
	}
	Anasazi::EpetraMultiVec Av(M_matrix->OperatorRangeMap(), M_indexNortol), u(M_matrix->OperatorRangeMap(), M_indexNortol);
	Anasazi::EpetraMultiVec* evecs = dynamic_cast<Anasazi::EpetraMultiVec* >(sol.Evecs->CloneViewNonConst( index ));

	Teuchos::SerialDenseMatrix<int,double> S(M_indexNortol,M_indexNortol);
	A->Apply( *evecs, Av );
	Av.MvNorm( tempnrm );
	for (int i=0; i<M_indexNortol; i++)
	{
		S(i,i) = one/tempnrm[i];
	};
	u.MvTimesMatAddMv( one, Av, S, zero );

	//
	// Compute direct residuals : || Av - sigma*u ||
	//

	for (int i=0; i<M_indexNortol; i++)
	{
		S(i,i) = evals[i].realpart;
	}

	Av.MvTimesMatAddMv( -one, u, S, one );
	Av.MvNorm( directnrm );
	if (M_comm->MyPID()==0)
	{
		std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
		std::cout<<std::setw(16)<<"Singular Value"
				<<std::setw(20)<<"Direct Residual"
				<<std::endl;
		std::cout<<"------------------------------------------------------"<<std::endl;
		for (int i=0; i<numev; i++)
		{
			std::cout<<std::setw(16)<<evals[i].realpart
					<<std::setw(20)<<directnrm[i]
					                           <<std::endl;
		}
		std::cout<<"------------------------------------------------------"<<std::endl;
	}

	//initializeV(M_indexNortol);

	std::vector<GIDtype> column_indices;
	for ( GIDtype i = 0; i < (GIDtype)(M_indexNortol); ++i )
		column_indices.push_back(i);

	GIDtype* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	// ATTENZIONE!!!!!
	//	M_map_column_unique.reset ( new Epetra_Map( (GIDtype)(-1), static_cast<int> (column_indices.size() ), (GIDtype)(0), *M_comm) ); // prima
	M_map_column_unique.reset ( new Epetra_Map( (GIDtype)(-1), static_cast<int> (column_indices.size() ), pointerToDofs, (GIDtype)(0), *M_comm) ); // adesso

	matrix.reset ( new Epetra_FECrsMatrix (Copy, M_matrix->OperatorRangeMap(), M_indexNortol ) );

	int * FirstPointInElementList1 = NULL;

	FirstPointInElementList1 = u.Map().FirstPointInElementList();

	double **A_Pointers = u.Pointers();

	for ( GIDtype i = 0; i < (GIDtype)(u.Map().NumMyElements()); i++ )
	{
		for ( GIDtype ii = 0; ii < (GIDtype)(u.Map().ElementSize(i)); ii++ )
		{
			GIDtype iii;
			iii = (GIDtype)(FirstPointInElementList1[i])+ii;
			for ( GIDtype j = 0; j < M_indexNortol; j++ )
			{
				GIDtype irow = (GIDtype)(u.Map().GID(i));
				GIDtype icol = (GIDtype)(j);

				matrix->InsertGlobalValues ( 1, &irow, 1, &icol, &A_Pointers[j][iii] );
			}
		}
	}

	matrix->GlobalAssemble(*M_map_column_unique, M_matrix->OperatorRangeMap());

}
//=========================================================================================
//////////////////////////////// LAPACK VERSION ////////////////////////////////
template <typename GIDtype>
void
V_POD<GIDtype>::performPOD_vector ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const Real nortol )
{
	Timer myTimer(M_comm, 1.0, "s");
	myTimer.StartGlobalTimer();

	checkTolerance(nortol);

	// A = M_matrix^T * M_matrix; in futuro A = M_matrix^T * matrice_norma *  M_matrix;
	boost::shared_ptr<Epetra_FECrsMatrix> A (new Epetra_FECrsMatrix(Copy, M_matrix->OperatorDomainMap(), 100));
	A->Scale(0.0);
	if (M_useWeightMatrix)
	{
		boost::shared_ptr<Epetra_FECrsMatrix> A_tmp (new Epetra_FECrsMatrix(Copy, M_weightMatrix->OperatorRangeMap(), 100));
		EpetraExt::MatrixMatrix::Multiply (*M_weightMatrix, false, *M_matrix, false, *A_tmp, false);
		A_tmp->GlobalAssemble( M_matrix->OperatorDomainMap(), M_weightMatrix->OperatorRangeMap() );
		EpetraExt::MatrixMatrix::Multiply (*M_matrix, true, *A_tmp, false, *A, false);
	}
	else
	{
		EpetraExt::MatrixMatrix::Multiply (*M_matrix, true, *M_matrix, false, *A, false);
	}
	A->GlobalAssemble(M_matrix->OperatorDomainMap(),M_matrix->OperatorDomainMap());

	Teuchos::SerialDenseMatrix<int, double> eigenvectors( M_numParameters , M_numParameters );
	std::vector<double> eigenvalues(M_numParameters);

	// at this point A is distributed across the processor, we have to sum up its contributes
	// and send it to Proc 0
	GIDtype NumMyElements_target;
	if( M_comm->MyPID() == 0 )
		NumMyElements_target = M_numParameters;
	else
		NumMyElements_target = 0;

	Epetra_Map TargetMap((GIDtype)(-1),NumMyElements_target,(GIDtype)(0),*M_comm);
	Epetra_Export Exporter(A->OperatorDomainMap(),TargetMap);
	boost::shared_ptr<Epetra_FECrsMatrix> A_target ( new Epetra_FECrsMatrix ( Copy, TargetMap, M_numParameters ) );
	A_target->Export(*A,Exporter,Add);
	A_target->GlobalAssemble(TargetMap,TargetMap);

	// Only Proc 0 performs EIGENSOLVE
	int root = 0;

	if ( M_comm->MyPID() == root )
	{

		// convert A into a SerialDenseMatrix
		Teuchos::SerialDenseMatrix<int, double> A_dense( M_numParameters , M_numParameters );
		FECrs_to_SerialDense(A_target, A_dense );

        // run SVD on A_dense
		GESVD_wrapper(A_dense, eigenvectors, eigenvalues );
	}

	M_comm->Broadcast( eigenvectors.values(), eigenvectors.numRows()*eigenvectors.numCols(), root );
	M_comm->Broadcast( &eigenvalues[0], M_numParameters, root );

	// compute singular values
	//std::vector<double> singularvalues(M_numParameters);
	for (int i = 0; i < eigenvalues.size() ; i++ )
	{
		M_singularValues[i] = std::sqrt( eigenvalues[i] );
	}

	if ( M_comm->MyPID() == 0 )
	{
		std::cout << "\n==================================\n" << std::endl;
		std::cout << "singular values: \n" << std::endl;
		for (int i = 0; i < eigenvalues.size() ; i++ )
		{
			std::cout << "   " << M_singularValues[i]  << std::endl;
		}
	}

	// Check on nortol
    setIndexNortol( eigenvalues );

	// build column map for the basis V
	std::vector<GIDtype> column_indices;
	for ( GIDtype i = 0; i < (GIDtype)(M_indexNortol); ++i )
		column_indices.push_back(i);

	GIDtype* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	M_map_column_unique.reset ( new Epetra_Map( (GIDtype)(-1), static_cast<int> (column_indices.size() ), pointerToDofs, (GIDtype)(0), *M_comm) ); // adesso

	boost::shared_ptr<Epetra_FECrsMatrix> tmp (new Epetra_FECrsMatrix(Copy, M_matrix->OperatorDomainMap(), M_indexNortol));

	// fill matrix tmp
	for ( GIDtype i = 0; i < (GIDtype) (M_matrix->OperatorDomainMap().NumMyElements()); i++ )
	{
		for ( GIDtype j = 0; j < (GIDtype) (M_map_column_unique->NumMyElements()); j++ )
		{
			double value = eigenvectors( (int)(i),(int)(j) ) / M_singularValues[j];
			tmp->InsertGlobalValues ( 1, &i, 1, &j, &value );
		}
	}

	tmp->GlobalAssemble(*M_map_column_unique, M_matrix->OperatorDomainMap());
	//std::cout << *tmp;

	// fill basis V
	matrix.reset ( new Epetra_FECrsMatrix (Copy, M_matrix->OperatorRangeMap(), M_indexNortol ) );

	EpetraExt::MatrixMatrix::Multiply ( *M_matrix, false, *tmp, false, *matrix, false );

	matrix->GlobalAssemble(*M_map_column_unique, M_matrix->OperatorRangeMap());

    myTimer.StopGlobalTimer("performing POD");

}
//=========================================================================================
//////////////////////////////// LOBPCG VERSION ////////////////////////////////
template <typename GIDtype>
void
V_POD<GIDtype>::performPOD_matrix ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix, const Real nortol, const Real& LOBPCG_tol, const int& LOBPCG_maxRestarts  )
{
	checkTolerance(nortol);

	int blockSize = M_numParameters;
	int maxRestarts = LOBPCG_maxRestarts;
	int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
	Teuchos::LAPACK<int,double> lapack;
	double tol = LOBPCG_tol;
	std::string which = "LM";

	Teuchos::ParameterList MyPL;
	MyPL.set( "Verbosity", verbosity );
	MyPL.set( "Which", which );
	MyPL.set( "Block Size", blockSize );
	MyPL.set( "Num Blocks", M_numBlocks );
	MyPL.set( "Maximum Restarts", maxRestarts );
	MyPL.set( "Convergence Tolerance", tol );
	MyPL.set( "Full Ortho", false );
	MyPL.set( "Maximum Iterations", 1000 );


	typedef Anasazi::MultiVec<double> MV;
	typedef Anasazi::Operator<double> OP;

	//typedef Epetra_MultiVector MV; per matrici invece che operatori!
	//typedef Epetra_Operator OP;

	//std::cout << *M_matrix;


	// Create an Anasazi::EpetraMultiVec for an initial vector to start the solver.
	// Note: This needs to have the same number of columns as the blocksize.
	Teuchos::RCP<Anasazi::EpetraMultiVec> ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec(M_matrix->OperatorDomainMap(), blockSize) );
	ivec->MvRandom();

	Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix( *M_matrix ) );
	Teuchos::RCP<Anasazi::EpetraSymOp> Amat = Teuchos::rcp( new Anasazi::EpetraSymOp( A ) );

	/* Altra opzione
	boost::shared_ptr<Epetra_FECrsMatrix> C (new Epetra_FECrsMatrix(Copy, M_matrix->OperatorDomainMap(), 100));
	C->Scale(0.0);
	EpetraExt::MatrixMatrix::Multiply (*M_matrix, true, *M_matrix, false, *C, false);
	C->GlobalAssemble(M_matrix->OperatorDomainMap(),M_matrix->OperatorDomainMap());
	Teuchos::RCP<Epetra_FECrsMatrix> A = Teuchos::rcp( new Epetra_FECrsMatrix( *C ) );
	A->GlobalAssemble(M_matrix->OperatorDomainMap(),M_matrix->OperatorDomainMap());
	*/


	Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
			Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(Amat, ivec) );


	// Inform the eigenproblem that the matrix A is symmetric
	MyProblem->setHermitian(true);

	// Set the number of eigenvalues requested and the blocksize the solver should use
	MyProblem->setNEV( M_nev );

	// Inform the eigenproblem that you are finished passing it information
	bool boolret = MyProblem->setProblem();
	if (!boolret) {
		throw "Anasazi::BasicEigenproblem::setProblem() returned with error.";
	}

	// Initialize the Block Arnoldi solver
	Anasazi::LOBPCGSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

	// Solve the problem to the specified tolerances or length
	Anasazi::ReturnType returnCode = MySolverMgr.solve();
	if (returnCode != Anasazi::Converged && M_comm->MyPID() ==0) {
		std::cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << std::endl;
	}


	// Get the eigenvalues and eigenvectors from the eigenproblem
	Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
	std::vector<Anasazi::Value<double> > evals = sol.Evals;
	int numev = sol.numVecs;

	const double one = 1.0;
	const double zero = 0.0;

	//
	// Compute singular values which are the square root of the eigenvalues
	//

	if (M_comm->MyPID()==0)
	{
		std::cout<<"------------------------------------------------------"<<std::endl;
		std::cout<<"Computed Singular Values: "<<std::endl;
		std::cout<<"------------------------------------------------------"<<std::endl;
	}

	//
	// Check on nortol
	//

	Real denominator = 0.0;
	Real numerator = 0.0;

	for (int i=0; i<numev; i++)
	{
		denominator += evals[i].realpart;
	}


	if (M_nortol >= 1.0)
	{ M_indexNortol = M_nortol;
	}
	else
	{
		M_indexNortol = -1;
		int k = 0;
		Real check = 0.0;
		Real test = 1.0 - ( M_nortol*M_nortol );

		while ( M_indexNortol < 0 && k <= numev )
		{
			numerator = 0.0;

			for ( int i = 0; i < k; ++i)
			{
				numerator += evals[i].realpart;
			}

			check = numerator/denominator;

			if ( check >= test )
			{
				M_indexNortol = k;
			}
			else
			{
				++k;
			}
		}
	}

	for (int i=0; i<M_indexNortol; i++)
	{
		evals[i].realpart = Teuchos::ScalarTraits<double>::squareroot( evals[i].realpart );
	}
	//
	// Compute left singular vectors : u = Av/sigma
	//
	std::vector<double> tempnrm(M_indexNortol), directnrm(M_indexNortol);
	std::vector<int> index(M_indexNortol);
	for (int i=0; i<M_indexNortol; i++)
	{
		index[i] = i;
	}
	Anasazi::EpetraMultiVec Av(M_matrix->OperatorRangeMap(), M_indexNortol), u(M_matrix->OperatorRangeMap(), M_indexNortol);
	Anasazi::EpetraMultiVec* evecs = dynamic_cast<Anasazi::EpetraMultiVec* >(sol.Evecs->CloneViewNonConst( index ));

	Teuchos::SerialDenseMatrix<int,double> S(M_indexNortol,M_indexNortol);
	A->Apply( *evecs, Av );
	Av.MvNorm( tempnrm );
	for (int i=0; i<M_indexNortol; i++)
	{
		S(i,i) = one/tempnrm[i];
	};
	u.MvTimesMatAddMv( one, Av, S, zero );

	//
	// Compute direct residuals : || Av - sigma*u ||
	//

	for (int i=0; i<M_indexNortol; i++)
	{
		S(i,i) = evals[i].realpart;
	}

	Av.MvTimesMatAddMv( -one, u, S, one );
	Av.MvNorm( directnrm );
	if (M_comm->MyPID()==0)
	{
		std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
		std::cout<<std::setw(16)<<"Singular Value"
				<<std::setw(20)<<"Direct Residual"
				<<std::endl;
		std::cout<<"------------------------------------------------------"<<std::endl;
		for (int i=0; i<numev; i++)
		{
			std::cout<<std::setw(16)<<evals[i].realpart
					<<std::setw(20)<<directnrm[i]
					                           <<std::endl;
		}
		std::cout<<"------------------------------------------------------"<<std::endl;
	}


	std::vector<GIDtype> column_indices;
	for ( GIDtype i = 0; i < (GIDtype)(M_indexNortol); ++i )
		column_indices.push_back(i);

	GIDtype* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	// ATTENZIONE!!!!!
//	M_map_column_unique.reset ( new Epetra_Map( (GIDtype)(-1), static_cast<int> (column_indices.size() ), (GIDtype)(0), *M_comm) ); // prima
	M_map_column_unique.reset ( new Epetra_Map( (GIDtype)(-1), static_cast<int> (column_indices.size() ), pointerToDofs, (GIDtype)(0), *M_comm) ); // adesso

	matrix.reset ( new Epetra_FECrsMatrix (Copy, M_matrix->OperatorRangeMap(), M_indexNortol ) );

	int * FirstPointInElementList1 = NULL;

	FirstPointInElementList1 = u.Map().FirstPointInElementList();

	double **A_Pointers = u.Pointers();

	for ( GIDtype i = 0; i < (GIDtype)(u.Map().NumMyElements()); i++ )
	{
		for ( GIDtype ii = 0; ii < (GIDtype)(u.Map().ElementSize(i)); ii++ )
		{
			GIDtype iii;

			iii = (GIDtype)(FirstPointInElementList1[i])+ii;

			for ( GIDtype j = 0; j < M_indexNortol; j++ )
			{


				GIDtype irow = (GIDtype)(u.Map().GID64(i));
				GIDtype icol = (GIDtype)(j);

				matrix->InsertGlobalValues ( 1, &irow, 1, &icol, &A_Pointers[j][iii] );
			}
		}
	}

	matrix->GlobalAssemble(*M_map_column_unique, M_matrix->OperatorRangeMap());

}
//=========================================================================================
template <typename GIDtype>
void
V_POD<GIDtype>::initializeV( const int num_cols)
{
	std::vector<GIDtype> column_indices;
	for ( GIDtype i = 0; i < (GIDtype)(num_cols); ++i )
		column_indices.push_back(i);

	GIDtype* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	M_map_column_repeated.reset ( new Epetra_Map( (GIDtype)(-1), static_cast<int> (column_indices.size() ), pointerToDofs, (GIDtype)(0), *M_comm) );

	M_map_column_unique.reset ( new Epetra_Map ( Epetra_Util::Create_OneToOne_Map ( *M_map_column_repeated, false ) ) );
}
//=========================================================================================
} // namespace LifeV

#endif //  V_POD_H_
