#include <taurus/Core/HyperReduction/HyperReductionVector.hpp>
#include <taurus/Core/HyperReduction/DEIM.hpp>
#include <taurus/Core/HyperReduction/ImportRhsSnapshots.hpp>
#include <taurus/Core/ReducedBasis/V_POD.hpp>
#include "EpetraExt_HDF5.h"
#include <taurus/Core/Utilities/Utilities.hpp>
#include "Epetra_SerialDenseSVD.h"
#include <taurus/Core/Utilities/Timer.hpp>


//#include "Epetra_Vector.h"

using namespace std;

namespace LifeV
{

//=========================================================================
// Constructor //
HyperReductionVector::HyperReductionVector (const GetPot& dataFile, const std::string& dataFile_section, const commPtr_Type& communicator, const mapPtr_Type& map_rows) :
	M_comm             ( communicator ),
	M_datafile         ( dataFile ),
	M_datafileSection  ( dataFile_section ),
	M_NumSnapshots	   ( 0 ),
	M_NumPODbases      ( 0 ),
	M_map_rows		   ( map_rows )
{
	std::string data_field = M_datafileSection+"/prefix";
	M_prefix = M_datafile(data_field.c_str(), "F");
}
//=========================================================================
void HyperReductionVector::perform( )
{

	EpetraExt::HDF5 HDF5_importer(*M_comm);

	// Read Snapshots
	HDF5_importer.Open( M_datafile ("importer/name_input_file_system", "SystemSnapshots")+".h5" );
	HDF5_importer.Read("info", "num_"+M_prefix+"Snapshots", M_NumSnapshots);
	HDF5_importer.Close();

	ImportRhsSnapshots RHSreader (M_map_rows, M_comm, M_NumSnapshots, M_datafile, M_datafileSection );

	RHSreader.readSnapshots( );

	RHSreader.getMatrix( M_Snapshots );

	// POD
	if ( M_comm->MyPID() == 0 )
		std::cout << "\n POD on " << M_datafileSection << " snapshots ...\n";

	V_POD<int> podRHS ( M_comm, M_Snapshots, M_NumSnapshots );

	std::string data_field = M_datafileSection+"/nortol";
	podRHS.performPOD_vector ( M_PODbasis, M_datafile(data_field.c_str(), 1e-2) );

	M_SingularValues = podRHS.getSingularValues();

	M_NumPODbases  = podRHS.getNumCols();

	M_Snapshots.reset();// clear memory

	// DEIM
 	DEIM<int> DEIM_rhs (M_comm, M_PODbasis, M_NumPODbases );

	DEIM_rhs.performDEIM( );

	M_DEIM_Phi = DEIM_rhs.getInterpolationMatrix();

	M_deim_GID = DEIM_rhs.getIndices();

}
//=========================================================================
void HyperReductionVector::write(DenseHDF5& HDF5dense_exporter )
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
		HDF5dense_exporter.WriteVectorInt(this_prefix, "indices", M_deim_GID);
		HDF5dense_exporter.WriteVectorInt(this_prefix, "dofs", M_deim_GID);
		//HDF5dense_exporter.Close();
	}

}
//=========================================================================
void HyperReductionVector::project(const boost::shared_ptr<Epetra_FECrsMatrix>& V, const int numColsV, EpetraExt::HDF5& HDF5_exporter)
{

	Timer myTimer(M_comm, 1.0, "s");

	RTRed_Utils getV_rhs_column;
	//boost::shared_ptr<Epetra_CrsMatrix> V_rhs_transposed;
	//getV_rhs_column.transpose( M_PODbasis, V_rhs_transposed);

	boost::shared_ptr<Epetra_MultiVector> M_PODbasis_multivec;

	getV_rhs_column.FECrsMatrixToMultiVector(M_PODbasis, M_PODbasis_multivec, M_NumPODbases);

	boost::shared_ptr<Epetra_Vector> FN( new Epetra_Vector (V->OperatorDomainMap()) );
	//boost::shared_ptr<Epetra_Vector> Fh( new Epetra_Vector (M_PODbasis->OperatorRangeMap()) );

	for ( int i = 1; i <= M_NumPODbases; ++i )
	{
		/*
		int numEntries = V_rhs_transposed->NumMyEntries(i-1);
		double* srcValues= new double[numEntries];
		int* colIndices = new int[numEntries];
		getV_rhs_column.extractRow(V_rhs_transposed, i-1, srcValues, colIndices);

		Fh->PutScalar(0.0);
		Fh->ReplaceGlobalValues(numEntries, srcValues, colIndices);
		V->Multiply(true, *Fh, *FN);
		delete [] srcValues;
		delete [] colIndices;
		*/

		myTimer.StartTimer();
		V->Multiply(true, *(*M_PODbasis_multivec)(i-1), *FN);
		double elaps1 = myTimer.StopTimer();

		// sum contributes distributed over all processors and send it to Proc 0
		// see http://trilinos.sandia.gov/Trilinos10.6Tutorial.pdf page 22
		int NumMyElements_target;
		if( M_comm->MyPID() == 0 )
			NumMyElements_target = numColsV;
		else
			NumMyElements_target = 0;

		Epetra_Map TargetMap(-1,NumMyElements_target,0,*M_comm);
		Epetra_Export Exporter(V->OperatorDomainMap(),TargetMap);
		boost::shared_ptr<Epetra_Vector> FN_target ( new Epetra_Vector (TargetMap) );
		FN_target->Export(*FN,Exporter,Add);

		myTimer.StartTimer();
		HDF5_exporter.Write(M_prefix+"N_"+EpetraExt::toString(i), *FN_target);
		double elaps2 = myTimer.StopTimer();

		if (M_comm->MyPID()==0)
		{
			std::streamsize ss = cout.precision();
			cout << std::setprecision(2) << std::scientific;
			cout << "Vector " << M_prefix+"N_" << i << " computed in " << elaps1
								   << " s and exported in " <<  elaps2 << " s\n";
			cout.unsetf(ios_base::fixed);
			cout.precision(ss);
		}
	}

	if ( M_comm->MyPID() == 0 )
	{
		cout << "============================\n";
	}

}
//=========================================================================
void HyperReductionVector::DEIMjacobianLeftProjection(const Epetra_CrsMatrix* V, boost::shared_ptr<Epetra_FECrsMatrix>& matrix)
{
    
    //buildMaskMatrix();
    
    // M_DEIM_Phi  = P^T * PHI
    
    // invert
    Epetra_SerialDenseSVD my_solver;
    
    my_solver.SetMatrix( M_DEIM_Phi );
    my_solver.Invert();
    
    Epetra_SerialDenseMatrix* DEIM_Phi_inv;
    
    DEIM_Phi_inv = my_solver.InvertedMatrix();
    
    boost::shared_ptr<Epetra_FECrsMatrix> DEIM_Phi_inv_Crs;
    
    SerialDense_to_FECrs( DEIM_Phi_inv, DEIM_Phi_inv_Crs );
    
    // matrix = V^T * DEIM_Phi_inv_Crs
    
    Epetra_FECrsMatrix *tmp;
    
    tmp = new Epetra_FECrsMatrix ( Copy, V->OperatorDomainMap(), M_NumPODbases );
    int errCode = EpetraExt::MatrixMatrix::Multiply ( *V, true, *M_PODbasis, false, *tmp, false );
    tmp->GlobalAssemble(M_PODbasis->OperatorDomainMap(),V->OperatorDomainMap());// attenzione!! check
    
    matrix.reset ( new Epetra_FECrsMatrix ( Copy, V->OperatorDomainMap(), M_NumPODbases ) );
    errCode = EpetraExt::MatrixMatrix::Multiply ( *tmp, false, *DEIM_Phi_inv_Crs, false, *matrix, false );
    matrix->GlobalAssemble(DEIM_Phi_inv_Crs->OperatorDomainMap(), V->OperatorDomainMap());// attenzione!! check
    
}
//=========================================================================
void HyperReductionVector::DEIMjacobianLeftProjection(const boost::shared_ptr<Epetra_FECrsMatrix>& V, EpetraExt::HDF5& HDF5_exporter, std::string dataset)
{

	Timer myTimer(M_comm, 1.0, "s");
	myTimer.StartTimer();

	boost::shared_ptr<Epetra_FECrsMatrix> matrix;

    // invert
    Epetra_SerialDenseSVD my_solver;

    my_solver.SetMatrix( M_DEIM_Phi );
    my_solver.Invert();

    Epetra_SerialDenseMatrix* DEIM_Phi_inv;

    DEIM_Phi_inv = my_solver.InvertedMatrix();

    //std::cout << "\n\n ==== DEIM_Phi_inv: \n" <<  *DEIM_Phi_inv << std::endl;

    boost::shared_ptr<Epetra_FECrsMatrix> DEIM_Phi_inv_Crs;

    SerialDense_to_FECrs( DEIM_Phi_inv, DEIM_Phi_inv_Crs );

    //std::cout << "\n\n ==== DEIM_Phi_inv_Crs: \n" <<  *DEIM_Phi_inv_Crs << std::endl;

    Epetra_FECrsMatrix *tmp;

    tmp = new Epetra_FECrsMatrix ( Copy, V->OperatorDomainMap(), M_NumPODbases );
    EpetraExt::MatrixMatrix::Multiply ( *V, true, *M_PODbasis, false, *tmp, false );
    tmp->GlobalAssemble(M_PODbasis->OperatorDomainMap(),V->OperatorDomainMap());

    matrix.reset ( new Epetra_FECrsMatrix ( Copy, V->OperatorDomainMap(), M_NumPODbases ) );
    EpetraExt::MatrixMatrix::Multiply ( *tmp, false, *DEIM_Phi_inv_Crs, false, *matrix, false );
    matrix->GlobalAssemble(M_PODbasis->OperatorDomainMap(), V->OperatorDomainMap());

    //std::cout << "\n\n ==== DEIMjacobianLeftProjection: \n" <<  *matrix << std::endl;

    // sum contributes distributed over all processors and send it to Proc 0
    // see http://trilinos.sandia.gov/Trilinos10.6Tutorial.pdf page 22
    int NumMyElements_Ctarget;
    if( M_comm->MyPID() == 0 )
    	NumMyElements_Ctarget = M_NumPODbases;
    else
    	NumMyElements_Ctarget = 0;

    int NumMyElements_Rtarget;
    if( M_comm->MyPID() == 0 )
    	NumMyElements_Rtarget = V->NumMyCols();
    else
    	NumMyElements_Rtarget = 0;

    Epetra_Map RTargetMap(-1,NumMyElements_Rtarget,0,*M_comm);
    Epetra_Map CTargetMap(-1,NumMyElements_Ctarget,0,*M_comm);

    Epetra_Export Exporter(matrix->OperatorRangeMap(),RTargetMap);
    boost::shared_ptr<Epetra_FECrsMatrix> matrix_target ( new Epetra_FECrsMatrix ( Copy, RTargetMap, M_NumPODbases ) );
    matrix_target->Export(*matrix,Exporter,Add);
    matrix_target->GlobalAssemble(CTargetMap,RTargetMap);

    HDF5_exporter.Write(dataset, *matrix_target);

    double elapsed_time = myTimer.StopTimer();
    if (M_comm->MyPID()==0)
    {
    	std::cout << "\nDEIMjacobianLeftProjection Matrix computed and exported in " << elapsed_time << "s\n\n";
    	std::cout << "============================\n";
    }

}
//=========================================================================
void HyperReductionVector::SerialDense_to_FECrs(const Epetra_SerialDenseMatrix* A_dense, boost::shared_ptr<Epetra_FECrsMatrix>& A )
{
    
    A.reset ( new Epetra_FECrsMatrix(Copy, M_PODbasis->OperatorDomainMap(), M_NumPODbases) );
    
    for ( int i = 0; i < A_dense->M(); i++ )
    {
        for ( int j = 0; j < A_dense->N(); j++ )
        {
            double value = (*A_dense)(i,j);
            //value  = value / M_comm->NumProc();
            A->InsertGlobalValues ( 1, &i, 1, &j, &value ) ;
        }
    }
    
    A->GlobalAssemble(M_PODbasis->OperatorDomainMap(), M_PODbasis->OperatorDomainMap());

}
//=========================================================================
void HyperReductionVector::buildMaskMatrix(boost::shared_ptr<Epetra_FECrsMatrix>& P)
{
 
    P.reset ( new Epetra_FECrsMatrix(Copy, M_PODbasis->OperatorRangeMap(), 1) );

    for (int j = 0 ; j < M_deim_GID.size(); ++j)
    {
        double value = 1.0;//
        value = 1.0 / M_comm->NumProc();
        P->InsertGlobalValues ( 1, &M_deim_GID[j], 1, &j, &value ) ;
    }
    P->GlobalAssemble(M_PODbasis->OperatorDomainMap(), M_PODbasis->OperatorRangeMap());

}
//=========================================================================
} // end namespace LifeV
