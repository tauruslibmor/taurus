#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <taurus/Core/Utilities/Timer.hpp>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

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

using namespace LifeV;
using namespace std;


void GESVD_wrapper(const Teuchos::SerialDenseMatrix<int,double>& A_dense, Teuchos::SerialDenseMatrix<int, double>& eigenvectors, std::vector<double>& eigenvalues )
{
	// for a description of GESVD inputs see
	// http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing.html#ga84fdf22a62b12ff364621e4713ce02f2
	// and
	// https://trilinos.org/docs/dev/packages/epetra/doc/html/classEpetra__LAPACK.html
	// See also https://trilinos.org/docs/dev/packages/teuchos/doc/html/LAPACK_2cxx_main_8cpp-example.html#_a0

	int N = A_dense.numRows();
	Epetra_LAPACK lapack;
	int info, lwork;
	double Vt[ 1 ];
	lwork = EPETRA_MAX( 3*N + N, 5*N );
	std::vector<double> work( lwork );

	char JOBU  = 'A';
	char JOBVT = 'N';

	std::cout << "\n GESVD running ...\n";

	lapack.GESVD( JOBU, JOBVT, N, N, A_dense.values(), A_dense.stride(), &eigenvalues[0], eigenvectors.values(), N, Vt, 1, &work[0], &lwork, &info );

	//std::cout << "\n==================================" << std::endl;

	//std::cout << "eigenvalues: \n" << std::endl;
	//for (int i = 0; i < eigenvalues.size() ; i++ )
	//{
	//	std::cout << "   " << eigenvalues[i]  << std::endl;
	//}

}

int
main ( int argc, char** argv )
{
    bool verbose (false);
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        verbose = true;
    }
#else
    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm () );
    verbose = true;
#endif

    int N = 10;
    bool parallel = true;

    if ( argc > 1 ) {
        N = atoi( argv[1] );
        parallel = atoi( argv[2] );
    }

    if (parallel)
    {
    	Teuchos::SerialDenseMatrix<int, double> A_dense(N,N);
    	A_dense.random();

    	Teuchos::SerialDenseMatrix<int, double> eigenvectors( N , N );
    	std::vector<double> eigenvalues(N);

    	Timer myTimer(Comm, 1.0, "s");
    	myTimer.StartTimer();
    	GESVD_wrapper(A_dense, eigenvectors, eigenvalues );
    	myTimer.StopTimer("solving Eigen-Problem with LAPACK GESVD");
    }
    else
    {
    	if (Comm->MyPID()==0)
    	{
    		Teuchos::SerialDenseMatrix<int, double> A_dense(N,N);
    		A_dense.random();

    		Teuchos::SerialDenseMatrix<int, double> eigenvectors( N , N );
    		std::vector<double> eigenvalues(N);

    		Timer myTimer(Comm, 1.0, "s");
    		myTimer.setSerial();
    		myTimer.StartTimer();
    		GESVD_wrapper(A_dense, eigenvectors, eigenvalues );
    		myTimer.StopTimer("solving Eigen-Problem with LAPACK GESVD");
    	}
    }

#ifdef HAVE_MPI
    if (verbose)
    {
        std::cout << "\nMPI Finalization\n" << std::endl;
    }
    MPI_Finalize();
#endif
    return ( EXIT_SUCCESS );
}
