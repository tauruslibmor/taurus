#include <taurus/CSM/Models/StructuralSolverTierIII.hpp>
#include <taurus/Core/HyperReduction/HyperReductionVector.hpp>
#include <taurus/Core/ReducedBasis/ReductionVector.hpp>
#include <taurus/Core/Utilities/Timer.hpp>

using namespace std;
using namespace std::chrono;

namespace LifeV
{

//=========================================================================
// Constructor //
StructuralSolverTierIII::StructuralSolverTierIII (const GetPot& dataFile, const commPtr_Type& communicator) :
			M_comm             ( communicator ),
			M_datafile         ( dataFile ),
			M_verbose		   ( false ),
			M_export_flag      ( true ),
			M_flagV			   ( false ),
			M_flagGuess 	   ( false ),
			M_pressureLoad	      (0.0),
			M_backtracingIter   ( 0 ),
			M_backtracingFactor ( 1.0 ),
			M_num_basis          ( 0 ),
			M_residual_Norm0     (0.0)
{
	if ( M_comm->MyPID() == 0 )
	{
		M_verbose = true;
	}
	M_Timer.reset ( new Timer (M_comm, 1e+3, "ms") );
	M_Timer->setPrecision(1);
	M_HDF5_importer.reset (new EpetraExt::HDF5(*M_comm) );
}
//=========================================================================
void StructuralSolverTierIII::SetVerbosity(const int& verbosity)
{
	switch (verbosity)
	{
	case (0):
		M_verbose = false;
		M_Timer->verbose(false);
		break;

	case (1):
		M_Timer->verbose(false);
	    break;

	default:
		M_Timer->verbose(true);
	}
}
//=========================================================================
void StructuralSolverTierIII::Setup(const meshPtr_Type& localMeshPtr, bcPtr_Type bc, bcPtr_Type bcMark, bool export_flag)
{
	M_Timer->StartTimer();
	M_bc		     = bc;
	M_bc_markingDOFs = bcMark;
	M_mesh = localMeshPtr;
	M_export_flag = export_flag;

	M_StructuralModel = M_datafile ("structure/model", "LinearElasticity");

	// SPAZI ELEMENTI FINITI
	M_FEspaceOrder = M_datafile ("finite_element/degree", "P1");

	M_dispFESpace.reset (new uSpaceStd_Type (M_mesh, M_FEspaceOrder, 3, M_comm) );
	M_FESpace_scalar.reset (new uSpaceStd_Type (M_mesh, M_FEspaceOrder, 1, M_comm) );
	M_dispETFESpace.reset ( new uSpaceETA_Type ( M_mesh, & ( M_dispFESpace->refFE() ), & ( M_dispFESpace->fe().geoMap() ), M_comm ) );

	M_solution.reset ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
	M_solution->zero();

	// BCHandler to mark dofs
	M_bc_markingDOFs->bcUpdate( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );
	M_bc->bcUpdate( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );

	// Marking dofs
	M_extractor.reset ( new DOF_Extractor ( M_dispFESpace ) ) ;
	M_extractor->setup( M_bc_markingDOFs );

	M_JacobianApproximation    = M_datafile ("hyperreduction/jacobian_approximation", "DEIM");
	M_SplitInternalForces      = M_datafile ("hyperreduction/split_internal_force", false);

	if ( M_verbose )
	{
		std::cout << "\n\n*************** CSM TierIII Solver  *********************";
		std::cout << "\n** Jacobian Approximation Strategy  = " << M_JacobianApproximation;
		std::cout << "\n** Split Internal Forces            = " << std::boolalpha << M_SplitInternalForces;
		std::cout << "\n******************************************************\n\n";
	}

	// Solid Model
	M_Solid_model.reset ( StructuralModel::StructuralModelFactory::instance().createObject (M_datafile("structure/model","none")));
	M_Solid_model->setFESpaces ( M_dispFESpace, M_dispETFESpace );
	M_Solid_model->setUp(M_datafile);

	M_pressureLoad = M_datafile ("structure/pressure_load", 0.0);

	M_backtracingIter    = M_datafile("newton/online_backtracing_iterations", 0);
	M_backtracingFactor  = M_datafile("newton/online_backtracing_factor", 1.0);

	M_Timer->StopTimer("setup CSM TierIII solver");

	load_V();
	load_HyperData();
}
//=========================================================================
void StructuralSolverTierIII::Setup(const meshPtr_Type& localMeshPtr, bcPtr_Type bc, bcPtr_Type bcMark, Epetra_CrsMatrix* sol_basis, bool export_flag)
{
	M_V = sol_basis;
	M_flagV = true;
	Setup(localMeshPtr, bc, bcMark, export_flag);
}
//=========================================================================
void StructuralSolverTierIII::SetPressureLoad(const double load)
{
	M_pressureLoad = load;
}
//=========================================================================
void StructuralSolverTierIII::setInitialGuess(const vectorPtr_Type& guess)
{
	M_flagGuess = true;
	M_NewtonGuess.reset (new VectorEpetra ( M_dispFESpace->map(), Unique ) );

	*M_NewtonGuess = *guess;
}
//=========================================================================
double StructuralSolverTierIII::computeFEMResidualNorm()
{
	M_Timer->StartTimer();

	vectorPtr_Type residual( new vector_Type ( M_dispFESpace->map(), Unique ) );

	M_Solid_model->evaluate_residual(M_solution, residual);

	bcManageRhs ( *residual, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );

	M_Timer->StopTimer("computign FEM residual norm");

	return residual->norm2();
}
//=========================================================================
double StructuralSolverTierIII::computeFEMResidualNorm(const vectorPtr_Type& sol)
{
	M_Timer->StartTimer();

	vectorPtr_Type residual( new vector_Type ( M_dispFESpace->map(), Unique ) );

	M_Solid_model->evaluate_residual(sol, residual);

	bcManageRhs ( *residual, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );

	M_Timer->StopTimer("computing FEM residual norm");

	return residual->norm2();
}
//=========================================================================
void StructuralSolverTierIII::load_V()
{
	// Load Matrix V (if it's the case)

	// open ROM HDF5 file manager
	M_HDF5_importer->Open( M_datafile ("importer/file_name", "ROM.h5") );

	M_HDF5_importer->Read("info", "num_basis", M_num_basis);

	std::vector<int> column_indices;
	for ( int i = 0; i < M_num_basis; ++i )
		column_indices.push_back(i);

	int* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	M_map_column_V.reset ( new MapEpetra ( -1, static_cast<int> (column_indices.size() ), pointerToDofs, M_comm ) );
	M_map_col.reset(new Epetra_Map(*M_map_column_V->map(Unique)) );

	if (!M_flagV)
	{
		if (M_verbose)
			cout << "\n Loading matrix V ... \n";

		M_Timer->StartTimer();
		M_HDF5_importer->Read ( "V", *M_map_col, *M_extractor->getMapUnmarked()->map(Unique), M_V );
		M_Timer->StopTimer("loading matrix V");
	}
}
//=========================================================================
void StructuralSolverTierIII::load_HyperData()
{

	M_Timer->StartTimer();

	M_h5_resIiso_prefix    = "/"+M_datafile("residual_internal_force_iso/prefix", "F")+"/";
	M_h5_resIvol_prefix    = "/"+M_datafile("residual_internal_force_vol/prefix", "F")+"/";
	M_h5_resE_prefix       = "/"+M_datafile("residual_external_force/prefix", "F")+"/";

	M_HDF5dense_importer.Open( M_datafile("rom_exporter/hyper_output_file", "ROM_Hyper")+".h5" );

	M_HDF5dense_importer.ReadVectorInt("/ReducedMesh/", "dofs", M_deim_dofs_GID);

	M_HDF5dense_importer.ReadVectorInt("/ReducedMesh/", "dofs_attached", M_deim_dofs_attached);

	M_HDF5dense_importer.ReadEpetraSerialDenseMatrix(M_h5_resIiso_prefix, "PHI", M_PHI_rhsIiso);

	M_HDF5dense_importer.ReadEpetraSerialDenseMatrix(M_h5_resIvol_prefix, "PHI", M_PHI_rhsIvol);

	M_HDF5dense_importer.ReadEpetraSerialDenseMatrix(M_h5_resE_prefix, "PHI", M_PHI_rhsE);

	M_HDF5dense_importer.ReadVectorInt(M_h5_resIiso_prefix, "indices", M_deim_rhsIiso_GID);

	M_HDF5dense_importer.ReadVectorInt(M_h5_resIvol_prefix, "indices", M_deim_rhsIvol_GID);

	M_HDF5dense_importer.ReadVectorInt(M_h5_resE_prefix, "indices", M_deim_rhsE_GID);

	// --------------------------------------------------------
	// iso forces
	int numColumnsPHI_Iiso = M_deim_rhsIiso_GID.size();
	std::vector<int> column_indicesResIso;
	for ( int i = 0; i < numColumnsPHI_Iiso; ++i )
		column_indicesResIso.push_back(i);

	int* pointerToDofsRes_Iso (0);
	if (column_indicesResIso.size() > 0)
		pointerToDofsRes_Iso = &column_indicesResIso[0];

	boost::shared_ptr<MapEpetra> map_column_LV_iso ( new MapEpetra ( -1, static_cast<int> (column_indicesResIso.size() ), pointerToDofsRes_Iso, M_comm ) );

	M_map_col_LV_iso.reset (new Epetra_Map (*map_column_LV_iso->map(Unique)));

	M_HDF5_importer->Read("DEIMjacobianLeftProjection_iso", *M_map_col_LV_iso, *M_map_col, M_LV_iso );

	// build mask matrix
	P_DEIM_iso.reset ( new Epetra_FECrsMatrix(Copy, M_V->OperatorRangeMap(), 1) );

	for (int j = 0 ; j < numColumnsPHI_Iiso; ++j)
	{
		double value = 1.0;//
		value = 1.0 / M_comm->NumProc();
		P_DEIM_iso->InsertGlobalValues ( 1, &M_deim_rhsIiso_GID[j], 1, &j, &value ) ;
	}
	P_DEIM_iso->GlobalAssemble(M_LV_iso->OperatorDomainMap(), M_V->OperatorRangeMap());

	// load resI vectors from hdf5
	Epetra_MultiVector* importedVector;
	for ( int i = 1; i <= numColumnsPHI_Iiso; ++i )
	{
		M_HDF5_importer->Read("resI_isoN_"+EpetraExt::toString(i), *M_map_col, importedVector, false );
		M_resIiso_n_imported.push_back(*importedVector);
	}

	// --------------------------------------------------------
	// vol forces
	int numColumnsPHI_Ivol = M_deim_rhsIvol_GID.size();
	std::vector<int> column_indicesResVol;
	for ( int i = 0; i < numColumnsPHI_Ivol; ++i )
		column_indicesResVol.push_back(i);

	int* pointerToDofsRes_Vol (0);
	if (column_indicesResVol.size() > 0)
		pointerToDofsRes_Vol = &column_indicesResVol[0];

	boost::shared_ptr<MapEpetra> map_column_LV_vol ( new MapEpetra ( -1, static_cast<int> (column_indicesResVol.size() ), pointerToDofsRes_Vol, M_comm ) );

	M_map_col_LV_vol.reset (new Epetra_Map (*map_column_LV_vol->map(Unique)));

	M_HDF5_importer->Read("DEIMjacobianLeftProjection_vol", *M_map_col_LV_vol, *M_map_col, M_LV_vol );

	// build mask matrix
	P_DEIM_vol.reset ( new Epetra_FECrsMatrix(Copy, M_V->OperatorRangeMap(), 1) );

	for (int j = 0 ; j < numColumnsPHI_Ivol; ++j)
	{
		double value = 1.0;//
		value = 1.0 / M_comm->NumProc();
		P_DEIM_vol->InsertGlobalValues ( 1, &M_deim_rhsIvol_GID[j], 1, &j, &value ) ;
	}
	P_DEIM_vol->GlobalAssemble(M_LV_vol->OperatorDomainMap(), M_V->OperatorRangeMap());

	// load resI vectors from hdf5
	for ( int i = 1; i <= numColumnsPHI_Ivol; ++i )
	{
		M_HDF5_importer->Read("resI_volN_"+EpetraExt::toString(i), *M_map_col, importedVector, false );
		M_resIvol_n_imported.push_back(*importedVector);
	}

	// --------------------------------------------------------
	// external forces
	for ( int i = 1; i <= M_deim_rhsE_GID.size(); ++i )
	{
		M_HDF5_importer->Read("resEN_"+EpetraExt::toString(i), *M_map_col, importedVector, false );
		M_resE_n_imported.push_back( *importedVector );
	}

	// --------------------------------------------------------

	// Close HDF5 file
	M_HDF5dense_importer.Close();
	M_Timer->StopTimer("dense HDF5");
}
//=========================================================================
void StructuralSolverTierIII::Solve(const std::vector<double>& parameters, meshSubPtr_Type& meshSub, int configuration_number)
{

	if (M_JacobianApproximation=="MDEIM")
		SolveMDEIM(parameters, meshSub, configuration_number);

	if (M_JacobianApproximation=="DEIM")
	{
		if (M_SplitInternalForces)
			SolveDEIM_split(parameters, meshSub, configuration_number);
		else
			SolveDEIM(parameters, meshSub, configuration_number);
	}

}
//=========================================================================
void StructuralSolverTierIII::updateExternalForces()
{
	M_Timer->StartTimer();

	vectorPtr_Type residual_ext( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type rhs_ext_I;

	M_Solid_model->evaluate_external_forces(M_solution, residual_ext, M_meshSub);
	M_extractor->restrictVector( residual_ext, "unmarked", rhs_ext_I );

	DEIM<int> DEIM_rhsE(M_comm, M_PHI_rhsE, M_deim_rhsE_GID);
	Epetra_SerialDenseVector rhs_DEIM_E;
	DEIM_rhsE.extractEntriesVector(rhs_ext_I, rhs_DEIM_E);
	Epetra_SerialDenseVector theta_vecE;
	DEIM_rhsE.solveOnline(rhs_DEIM_E, theta_vecE);

	int size_theta_vecE = theta_vecE.Length();

	M_resE_n.reset ( new Epetra_MultiVector(*M_map_col, 1) );
	M_resE_n->PutScalar(0.0);

	for ( int i = 1; i <= size_theta_vecE; ++i )
	{
		M_resE_n->Update( theta_vecE[i-1], M_resE_n_imported[i-1], 1.0);
	}

	M_Timer->StopTimer("assembling external forces");

}
//=========================================================================
void StructuralSolverTierIII::compute_residual(const vectorPtr_Type& x, multivecPtr_Type& residual)
{
	vectorPtr_Type residual_int( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_int_iso( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_int_vol( new vector_Type ( M_dispFESpace->map(), Unique ) );

	vectorPtr_Type rhs_I;
	vectorPtr_Type rhs_iso_I;
	vectorPtr_Type rhs_vol_I;
	vectorPtr_Type rhs_int_I;

	// FEM assembly on reduced mesh
	M_Timer->StartTimer();
	M_Solid_model->evaluate_internal_forces_iso(x, residual_int_iso, M_meshSub);
	M_Solid_model->evaluate_internal_forces_vol(x, residual_int_vol, M_meshSub);
	M_Timer->StopTimer("assembling internal forces");

	*residual_int_vol *= -1.0;
	*residual_int_iso *= -1.0;

	// extract internal dofs
	M_Timer->StartTimer();
	M_extractor->restrictVectorFast( residual_int_iso, M_deim_dofs_attached, rhs_iso_I );
	M_extractor->restrictVectorFast( residual_int_vol, M_deim_dofs_attached, rhs_vol_I );
	M_Timer->StopTimer("restricting residual to internal nodes");

	// solve Interpolation problem Ax = b for iso contribute
	M_Timer->StartTimer();
	DEIM<int> DEIM_rhsI_iso(M_comm, M_PHI_rhsIiso, M_deim_rhsIiso_GID);
	Epetra_SerialDenseVector rhs_DEIM_Iiso;
	DEIM_rhsI_iso.extractEntriesVector(rhs_iso_I, rhs_DEIM_Iiso);

	Epetra_SerialDenseVector theta_vecIiso;
	DEIM_rhsI_iso.solveOnline(rhs_DEIM_Iiso, theta_vecIiso);

	M_Timer->StopTimer("DEIM iso rhs");

	// solve Interpolation problem Ax = b for vol contribute
	M_Timer->StartTimer();
	DEIM<int> DEIM_rhsI_vol(M_comm, M_PHI_rhsIvol, M_deim_rhsIvol_GID);
	Epetra_SerialDenseVector rhs_DEIM_Ivol;
	DEIM_rhsI_vol.extractEntriesVector(rhs_vol_I, rhs_DEIM_Ivol);

	Epetra_SerialDenseVector theta_vecIvol;
	DEIM_rhsI_vol.solveOnline(rhs_DEIM_Ivol, theta_vecIvol);

	M_Timer->StopTimer("DEIM vol rhs");

	// sum affine contributes of iso and vol contributes
	M_Timer->StartTimer();
	// iso
	int size_theta_vecIiso = theta_vecIiso.Length();
	M_resIiso_n.reset ( new Epetra_MultiVector(*M_map_col, 1) );
	M_resIiso_n->PutScalar(0.0);

	for ( int i = 1; i <= size_theta_vecIiso; ++i )
	{
		M_resIiso_n->Update( theta_vecIiso[i-1], M_resIiso_n_imported[i-1], 1.0);
	}

	// vol
	int size_theta_vecIvol = theta_vecIvol.Length();
	M_resIvol_n.reset ( new Epetra_MultiVector(*M_map_col, 1) );
	M_resIvol_n->PutScalar(0.0);

	for ( int i = 1; i <= size_theta_vecIvol; ++i )
	{
		M_resIvol_n->Update( theta_vecIvol[i-1], M_resIvol_n_imported[i-1], 1.0);
	}

	// sum iso and vol final contributes
	residual.reset ( new Epetra_MultiVector ( *M_map_col, 1) );
	residual->Update( 1.0,  *M_resIiso_n,  1.0, *M_resIvol_n, 0.0);
	residual->Update( -1.0, *M_resE_n, 1.0);

	M_Timer->StopTimer("forming rhs reduced");
}
//=========================================================================
void StructuralSolverTierIII::jacobian_projection(const matrixPtr_Type& jac, const boost::shared_ptr<Epetra_FECrsMatrix>& mask, const Epetra_CrsMatrix* Left, Epetra_FECrsMatrix*& result)
{
	// jac = left * ( ( mask * jac ) * M_V )

	boost::shared_ptr<Epetra_FECrsMatrix> J_DEIM ( new Epetra_FECrsMatrix ( Copy, mask->OperatorDomainMap(), 50 ) );
	EpetraExt::MatrixMatrix::Multiply ( *mask, true, *jac->matrixPtr(), false, *J_DEIM, false );
	J_DEIM->GlobalAssemble( mask->OperatorRangeMap(), mask->OperatorDomainMap() );

	boost::shared_ptr<Epetra_FECrsMatrix> J_DEIM_V ( new Epetra_FECrsMatrix ( Copy, mask->OperatorDomainMap(), M_num_basis ) );
	EpetraExt::MatrixMatrix::Multiply ( *J_DEIM, false, *M_V, false, *J_DEIM_V, false );
	J_DEIM_V->GlobalAssemble( M_V->OperatorDomainMap(), mask->OperatorDomainMap() );

	result =  new Epetra_FECrsMatrix ( Copy, M_V->OperatorDomainMap(), M_num_basis ) ;
	EpetraExt::MatrixMatrix::Multiply ( *Left, false, *J_DEIM_V, false, *result, false );
	result->GlobalAssemble(M_V->OperatorDomainMap(),M_V->OperatorDomainMap());

}
//=========================================================================
void StructuralSolverTierIII::compute_jacobian_DEIM(const vectorPtr_Type& x, crsPtr_Type& jacobian)
{
	// FEM assembly on reduced mesh
	M_Timer->StartTimer();

	matrixPtr_Type jacobian_iso ( new matrix_Type ( M_dispETFESpace->map(), 100 ) );
	matrixPtr_Type jacobian_vol ( new matrix_Type ( M_dispETFESpace->map(), 100 ) );

	M_Solid_model->evaluate_jacobian_iso(x, jacobian_iso, M_meshSub);
	M_Solid_model->evaluate_jacobian_vol(x, jacobian_vol, M_meshSub);

	M_Timer->StopTimer("assembling jacobian on reduced mesh");

	// extract internal dofs
	M_Timer->StartTimer();

	matrixPtr_Type matrix_I_I_iso;
	matrixPtr_Type matrix_I_I_vol;

	M_extractor->restrictMatrixFast( jacobian_iso, M_deim_dofs_attached, matrix_I_I_iso );
	M_extractor->restrictMatrixFast( jacobian_vol, M_deim_dofs_attached, matrix_I_I_vol );

	M_Timer->StopTimer("restricting jacobian to internal nodes");

	// DEIM-Reduced Jacobian
	M_Timer->StartTimer();

	// iso contribute
	Epetra_FECrsMatrix *An_DEIM_iso;
	jacobian_projection(matrix_I_I_iso, P_DEIM_iso, M_LV_iso, An_DEIM_iso);

	// vol contribute
	Epetra_FECrsMatrix *An_DEIM_vol;
	jacobian_projection(matrix_I_I_vol, P_DEIM_vol, M_LV_vol, An_DEIM_vol);

	// sum contributes
	jacobian.reset ( new Epetra_CrsMatrix ( Copy, *M_map_col, M_num_basis, false ) );
	EpetraExt::MatrixMatrix::Add ( *An_DEIM_iso, false, 1.0, *jacobian, 1.0 );
	EpetraExt::MatrixMatrix::Add ( *An_DEIM_vol, false, 1.0, *jacobian, 1.0 );

	jacobian->FillComplete();

	M_Timer->StopTimer("forming DEIM-approximated jacobian");

}
//=========================================================================
void StructuralSolverTierIII::backtracing(const vectorPtr_Type& x_k, const vectorPtr_Type& delta, double norm_old, double& alpha)
{
	vectorPtr_Type x_kp1 ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
	x_kp1->zero();
	multivecPtr_Type fn;

	alpha = 1;
	double residual_Norm;
	int backtrack_iter = 0;
    
    *x_kp1 = *x_k + alpha * (*delta);
    compute_residual(x_kp1, fn);
    fn->Norm2(&residual_Norm);

	while (residual_Norm > norm_old && backtrack_iter < M_backtracingIter )
	{
		backtrack_iter = backtrack_iter + 1;

		alpha = alpha * M_backtracingFactor;

		*x_kp1 = *x_k + alpha * (*delta);

		compute_residual(x_kp1, fn);
		fn->Norm2(&residual_Norm);

		if (M_verbose)
			cout << "\n      Backtracking: Alpha = " << alpha << ", Relative residual norm = " << residual_Norm / M_residual_Norm0 << "\n";
	}
}
//=========================================================================
void StructuralSolverTierIII::SolveDEIM_split(const std::vector<double>& parameters, meshSubPtr_Type& meshSub, int configuration_number)
{
	StartGlobalTimer();

	M_meshSub = meshSub;// check

	// VETTORE SOLUZIONE
	vectorPtr_Type solution_rb ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
	if (M_flagGuess)
		*solution_rb = *M_NewtonGuess;
	else
		solution_rb->zero();

	// Set Exporter
	ExporterHDF5< mesh_Type > exporter(M_datafile, "exporter");
	if (M_export_flag)
	{
		exporter.setMeshProcId( M_mesh, M_comm->MyPID() );
		string FileName =  "SolutionTierIII"+EpetraExt::toString(configuration_number);
		exporter.setPrefix( FileName);
		exporter.setPostDir( "./" );
		exporter.addVariable ( ExporterData< mesh_Type >::VectorField, "SolutionTierIII", M_dispFESpace, solution_rb, UInt ( 0 ) );
	}

	// Pseudo time step for the exporter
	Real istance = 0.0;

	vectorPtr_Type sol ( new vector_Type (*M_map_column_V, Unique ) );
	vectorPtr_Type increment ( new vector_Type (*M_map_column_V, Unique ) );

	vectorPtr_Type increment_h ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
	increment_h->zero();

	// Material model
	if (M_export_flag)
		M_Solid_model->setExportStresses(exporter, M_FESpace_scalar);

	M_Solid_model->updateLoad(M_pressureLoad);

	// Assemble External forces once and for all
	updateExternalForces();

	//=========================================================================================
	// START NEWTON ITERATIONS
	int Newton_MaxIt = M_datafile("newton/maxiter",5);
	int Newton_It = 0;
	double Newton_tolerance = M_datafile("newton/tolerance",1e-6);
	double increment_Norm = Newton_tolerance + 1;
	double residual_Norm = 0;

	multivecPtr_Type fn;
	crsPtr_Type An_DEIM;

	M_bc->bcUpdate( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );
	bcManageRhs ( *solution_rb, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );

	compute_residual(solution_rb, fn);
	fn->Norm2(&M_residual_Norm0);

	if ( M_comm->MyPID() == 0 )
		std::cout << "\n\n========= Start Newton Iterations ==========\n\n";

	// Newton loop
	while (increment_Norm > Newton_tolerance && Newton_It < Newton_MaxIt)
	{
		Newton_It  = Newton_It + 1;
		increment->zero();
		increment_h->zero();
		sol->zero();

		vectorPtr_Type increment_h_I ( new VectorEpetra ( M_extractor->getMapUnmarked() ) );
		increment_h_I->zero();

		compute_jacobian_DEIM(solution_rb, An_DEIM);

		//------------------------------------------------------------
		// Solve RB problem
		M_Timer->StartTimer();

		Amesos_BaseSolver * Solver;
		Amesos Factory;
		std::string SolverType = M_datafile ("solver/method", "Klu");

		Epetra_LinearProblem Problem;
		// solve
		Problem.SetOperator(&(*An_DEIM));
		Problem.SetRHS(&(*fn));
		Problem.SetLHS(&sol->epetraVector());

		Solver = Factory.Create(SolverType, Problem);
		Solver->SymbolicFactorization();
		Solver->NumericFactorization();
		Solver->Solve();

		M_Timer->StopTimer("solving RB system");
		//------------------------------------------------------------

		M_V->Multiply (false, sol->epetraVector(), increment_h_I->epetraVector());
		M_extractor->sumUnmarkedIntoGlobal( increment_h_I, increment_h );

		double alpha = 1;
		backtracing(solution_rb, increment_h, residual_Norm, alpha);

		*increment_h = alpha * (*increment_h);
		*solution_rb +=  (*increment_h);

		increment_Norm = (increment_h->norm2()) / (solution_rb->norm2());

		compute_residual(solution_rb, fn);
		fn->Norm2(&residual_Norm);

		if ( M_verbose )
		{
			//std::cout << "\n**** Iteration " << Newton_It << ", Norm Increment = " <<  increment_Norm << ", Residual Norm = " << residual_Norm/M_residual_Norm0 << "\n";
			printf("\n**** Iteration %d, Norm Increment = %1.2e, Rel. Residual Norm = %1.2e\n", Newton_It, increment_Norm, residual_Norm/M_residual_Norm0);
			std::cout << "===================================================\n";
		}

	}// END NEWTON ITERATIONS

	//=========================================================================================
	// Post-processing
	M_Timer->StartTimer();

	if (M_export_flag)
	{
		M_Solid_model->computeElementStresses(*solution_rb);
		exporter.postProcess( istance );
		exporter.closeFile();
	}
	// Close HDF5 Files
	M_HDF5_importer->Close();

	M_Timer->StopTimer("PostProcessing");
	*M_solution = *solution_rb;
	StopGlobalTimer();

}
//=========================================================================
void StructuralSolverTierIII::SolveDEIM(const std::vector<double>& parameters, meshSubPtr_Type& meshSub, int configuration_number)
{
	StartGlobalTimer();
	//=========================================================================================
	// Hyper-reduction HDF5 file manager
	StartTimer();

	std::string h5_resI_prefix    = "/"+M_datafile("residual_internal_force/prefix", "F")+"/";
	std::string h5_resE_prefix    = "/"+M_datafile("residual_external_force/prefix", "F")+"/";

	DenseHDF5 HDF5dense_importer;

	HDF5dense_importer.Open( M_datafile("rom_exporter/hyper_output_file", "ROM_Hyper")+".h5" );

	std::vector<int> deim_dofs_GID;
	HDF5dense_importer.ReadVectorInt("/ReducedMesh/", "dofs", deim_dofs_GID);

	Epetra_SerialDenseMatrix PHI_rhsI;
	HDF5dense_importer.ReadEpetraSerialDenseMatrix(h5_resI_prefix, "PHI", PHI_rhsI);

	Epetra_SerialDenseMatrix PHI_rhsE;
	HDF5dense_importer.ReadEpetraSerialDenseMatrix(h5_resE_prefix, "PHI", PHI_rhsE);

	std::vector<int> deim_rhsI_GID;
	HDF5dense_importer.ReadVectorInt(h5_resI_prefix, "indices", deim_rhsI_GID);

	std::vector<int> deim_rhsE_GID;
	HDF5dense_importer.ReadVectorInt(h5_resE_prefix, "indices", deim_rhsE_GID);

	// Close HDF5 file
	HDF5dense_importer.Close();

	StopTimer("dense HDF5");

	//=========================================================================================
	mapPtr_Type map_rows;
	map_rows.reset ( new MapEpetra ( *M_extractor->getMapUnmarked() ) );

	// VETTORE SOLUZIONE
	vectorPtr_Type solution_rb ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
	solution_rb->zero();

	vectorPtr_Type increment_h ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
	increment_h->zero();

	// Set Exporter
	ExporterHDF5< mesh_Type > exporter(M_datafile, "exporter");
	if (M_export_flag)
	{
		exporter.setMeshProcId( M_mesh, M_comm->MyPID() );
		string FileName =  "SolutionTierIII"+EpetraExt::toString(configuration_number);
		exporter.setPrefix( FileName);
		exporter.setPostDir( "./" );
		exporter.addVariable ( ExporterData< mesh_Type >::VectorField, "SolutionTierIII", M_dispFESpace, solution_rb, UInt ( 0 ) );
	}

	matrixPtr_Type matrix_I_I;
	matrixPtr_Type matrix_I_D;
	vectorPtr_Type rhs_I;
	vectorPtr_Type rhs_D;
	vectorPtr_Type rhs_ext_I;
	vectorPtr_Type rhs_int_I;

	// Pseudo time step for the exporter
	Real istance = 0.0;

	//=============
	// Load Matrix V
	StartTimer();

	// ROM HDF5 file manager
	EpetraExt::HDF5 HDF5_importer(*M_comm);
	HDF5_importer.Open( M_datafile ("importer/file_name", "ROM.h5") );

	int num_basis;
	HDF5_importer.Read("info", "num_basis", num_basis);

	std::vector<int> column_indices;
	for ( int i = 0; i < num_basis; ++i )
		column_indices.push_back(i);

	int* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	boost::shared_ptr<MapEpetra> map_column_V ( new MapEpetra ( -1, static_cast<int> (column_indices.size() ), pointerToDofs, M_comm ) );

	Epetra_Map map_col(*map_column_V->map(Unique));

	Epetra_CrsMatrix* M_V;
	HDF5_importer.Read ( "V", map_col, *M_extractor->getMapUnmarked()->map(Unique), M_V );

	StopTimer("loading matrix V");

	vectorPtr_Type sol ( new vector_Type (*map_column_V, Unique ) );
	vectorPtr_Type increment ( new vector_Type (*map_column_V, Unique ) );
	//=============


	// Jacobian matrix and residual vector
	vectorPtr_Type residual( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_ext( new vector_Type ( M_dispFESpace->map(), Unique ) );

	// Material model
	if (M_export_flag)
		M_Solid_model->setExportStresses(exporter, M_FESpace_scalar);

	M_Solid_model->updateLoad(M_pressureLoad);

	//=========================================================================================
	// START NEWTON ITERATIONS
	int Newton_MaxIt = M_datafile("newton/maxiter",5);
	int Newton_It = 0;
	double Newton_tolerance = M_datafile("newton/tolerance",1e-6);
	double increment_Norm = Newton_tolerance + 1;
	double residual_Norm = 0;
	double residual_Norm0 = 0;
	double residual_Norm_old = 0;

	M_bc->bcUpdate( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );
	bcManageRhs ( *solution_rb, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );

	//=========================================================================================
	// Assemble External forces once and for all
	M_Solid_model->evaluate_external_forces(solution_rb, residual_ext, meshSub);
	M_extractor->restrictVector( residual_ext, "unmarked", rhs_ext_I );

	DEIM<int> DEIM_rhsE(M_comm, PHI_rhsE, deim_rhsE_GID);
	Epetra_SerialDenseVector rhs_DEIM_E;
	DEIM_rhsE.extractEntriesVector(rhs_ext_I, rhs_DEIM_E);
	Epetra_SerialDenseVector theta_vecE;
	DEIM_rhsE.solveOnline(rhs_DEIM_E, theta_vecE);

	int size_theta_vecE = theta_vecE.Length();
	Epetra_MultiVector resE_n ( *map_column_V->map(Unique), 1 );
	resE_n.PutScalar(0.0);

	Epetra_MultiVector* importedVector;
	for ( int i = 1; i <= size_theta_vecE; ++i )
	{
		HDF5_importer.Read("resEN_"+EpetraExt::toString(i), map_col, importedVector, false );
		resE_n.Update( theta_vecE[i-1], *importedVector, 1.0);
	}

	if ( M_comm->MyPID() == 0 )
		std::cout << "\n\n========= Start Newton Iterations ==========\n\n";


	//=========================================================================================
	// Structures for DEIM Jacobian approximation

	// column map
	int numColumnsPHI_I = deim_rhsI_GID.size();
	std::vector<int> column_indicesRes;
	for ( int i = 0; i < numColumnsPHI_I; ++i )
		column_indicesRes.push_back(i);

	int* pointerToDofsRes (0);
	if (column_indicesRes.size() > 0)
		pointerToDofsRes = &column_indicesRes[0];

	boost::shared_ptr<MapEpetra> map_column_LV ( new MapEpetra ( -1, static_cast<int> (column_indicesRes.size() ), pointerToDofsRes, M_comm ) );

	Epetra_Map map_col_LV(*map_column_LV->map(Unique));

	// load left projection matrix
	Epetra_CrsMatrix * LV;
	HDF5_importer.Read("DEIMjacobianLeftProjection", map_col_LV, map_col, LV );

	// build mask matrix
	boost::shared_ptr<Epetra_FECrsMatrix> P_DEIM;
	P_DEIM.reset ( new Epetra_FECrsMatrix(Copy, M_V->OperatorRangeMap(), 1) );

	for (int j = 0 ; j < numColumnsPHI_I; ++j)
	{
		double value = 1.0;//
		value = 1.0 / M_comm->NumProc();
		P_DEIM->InsertGlobalValues ( 1, &deim_rhsI_GID[j], 1, &j, &value ) ;
	}
	P_DEIM->GlobalAssemble(LV->OperatorDomainMap(), M_V->OperatorRangeMap());

	// load resI vectors from hdf5
	std::vector<Epetra_MultiVector> resI_n_imported; //( *map_column_V->map(Unique), 1 );
	for ( int i = 1; i <= numColumnsPHI_I; ++i )
	{
		HDF5_importer.Read("resIN_"+EpetraExt::toString(i), map_col, importedVector, false );
		resI_n_imported.push_back(*importedVector);
	}

	// Newton loop
	while (increment_Norm > Newton_tolerance && Newton_It < Newton_MaxIt)
	{
		matrixPtr_Type jacobian ( new matrix_Type ( M_dispETFESpace->map(), 100 ) );
		vectorPtr_Type residual_int( new vector_Type ( M_dispFESpace->map(), Unique ) );

		Newton_It  = Newton_It + 1;
		increment->zero();
		increment_h->zero();
		sol->zero();

		//=========================================================================================
		// ASSEMBLAGGIO SU MESH RIDOTTA
		StartTimer();
		M_Solid_model->evaluate_jacobian(solution_rb, jacobian, meshSub);
		//=========================================================================================
		// ORA DEVO ESTRARRE DALLA MATRICE LE COMPONENTI INTERNE,DI DIRICHLET E MISTE

		M_extractor->restrictMatrix( jacobian, "unmarked", "unmarked", matrix_I_I );
		M_extractor->restrictMatrix( jacobian, "unmarked", "marked",  matrix_I_D );
		StopTimer("assembling jacobian on reduced mesh");

		//=========================================================================================
		// ASSEMBLAGGIO RHS SU MESH RIDOTTA
		StartTimer();
		M_Solid_model->evaluate_internal_forces(solution_rb, residual_int, meshSub);
		bcManageRhs ( *residual_int, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );
		*residual_int *= -1.0;

		//=========================================================================================
		// ORA DEVO ESTRARRE DAL RHS LE COMPONENTI INTERNE E DI DIRICHLET

		M_extractor->restrictVector( residual_int, "unmarked", rhs_I );
		M_extractor->restrictVector( residual_int, "marked", rhs_D );

		vectorPtr_Type rhs_ridotto ( new VectorEpetra ( *rhs_I - (*matrix_I_D * (*rhs_D) ) ) );

		vectorPtr_Type increment_h_I ( new VectorEpetra ( rhs_I->map() ) );
		increment_h_I->zero();
		StopTimer("assembling internal forces");

		//=========================================================================================
		// solve Interpolation problem Ax = b for RHS_I

		StartTimer();
		DEIM<int> DEIM_rhsI(M_comm, PHI_rhsI, deim_rhsI_GID);
		Epetra_SerialDenseVector rhs_DEIM_I;
		DEIM_rhsI.extractEntriesVector(rhs_ridotto, rhs_DEIM_I);

		Epetra_SerialDenseVector theta_vecI;
		DEIM_rhsI.solveOnline(rhs_DEIM_I, theta_vecI);

		StopTimer("DEIM rhs");

		StartTimer();
		int size_theta_vecI = theta_vecI.Length();
		Epetra_MultiVector resI_n ( *map_column_V->map(Unique), 1 );
		resI_n.PutScalar(0.0);

		for ( int i = 1; i <= size_theta_vecI; ++i )
		{
			//HDF5_importer.Read("resIN_"+EpetraExt::toString(i), map_col, importedVector, false );
			resI_n.Update( theta_vecI[i-1], /**importedVector*/ resI_n_imported[i-1], 1.0);
		}

		Epetra_MultiVector fn ( *map_column_V->map(Unique), 1 );
		fn.Update( 1.0,  resI_n, -1.0, resE_n, 0.0);

		fn.Norm2(&residual_Norm);
		if (Newton_It == 1)
			residual_Norm0 = residual_Norm;

		StopTimer("forming rhs reduced");

		//=========================================================================================
		// DEIM-Reduced Jacobian
		StartTimer();

		boost::shared_ptr<Epetra_FECrsMatrix> J_DEIM ( new Epetra_FECrsMatrix ( Copy, P_DEIM->OperatorDomainMap(), 50 ) );
		int errCode = EpetraExt::MatrixMatrix::Multiply ( *P_DEIM, true, *matrix_I_I->matrixPtr(), false, *J_DEIM, false );
		J_DEIM->GlobalAssemble( P_DEIM->OperatorRangeMap(), P_DEIM->OperatorDomainMap() );

		boost::shared_ptr<Epetra_FECrsMatrix> J_DEIM_V ( new Epetra_FECrsMatrix ( Copy, P_DEIM->OperatorDomainMap(), num_basis ) );
		errCode = EpetraExt::MatrixMatrix::Multiply ( *J_DEIM, false, *M_V, false, *J_DEIM_V, false );
		J_DEIM_V->GlobalAssemble( M_V->OperatorDomainMap(), P_DEIM->OperatorDomainMap() );

		Epetra_FECrsMatrix *An_DEIM;

		An_DEIM =  new Epetra_FECrsMatrix ( Copy, M_V->OperatorDomainMap(), num_basis ) ;
		errCode = EpetraExt::MatrixMatrix::Multiply ( *LV, false, *J_DEIM_V, false, *An_DEIM, false );
		An_DEIM->GlobalAssemble(M_V->OperatorDomainMap(),M_V->OperatorDomainMap());

		StopTimer("forming DEIM-approximated jacobian");

		//=========================================================================================
		// Solve RB problem
		StartTimer();

		Amesos_BaseSolver * Solver;
		Amesos Factory;
		std::string SolverType = M_datafile ("solver/method", "Klu");

		Epetra_LinearProblem Problem;
		// solve
		Problem.SetOperator(An_DEIM);
		Problem.SetRHS(&fn);
		Problem.SetLHS(&sol->epetraVector());

		Solver = Factory.Create(SolverType, Problem);
		Solver->SymbolicFactorization();
		Solver->NumericFactorization();
		Solver->Solve();

		StopTimer("solving RB system");

		M_V->Multiply (false, sol->epetraVector(), increment_h_I->epetraVector());
		M_extractor->sumUnmarkedIntoGlobal( increment_h_I, increment_h );

		*solution_rb += *increment_h;

		matrix_I_I->zero();
		matrix_I_D->zero();
		rhs_I->zero();
		rhs_D->zero();
		rhs_ridotto->zero();

		increment_Norm = (increment_h->norm2()) / (solution_rb->norm2());


		if ( M_comm->MyPID() == 0 )
		{
			std::cout << "\n**** Iteration " << Newton_It << ", Norm Increment = " <<  increment_Norm << ", Residual Norm = " << residual_Norm/residual_Norm0 << "\n";
			std::cout << "===================================================\n";
		}

		residual_Norm_old = residual_Norm;

		//if (M_export_flag)
		//{
		//	M_Solid_model->computeElementStresses(*solution_rb);
		//	exporter.postProcess( istance );
		//	//exporter.closeFile();
		//	istance+=1;
		//}


	}// END NEWTON ITERATIONS


	//=========================================================================================
	//=========================================================================================
	// Post-processing
	StartTimer();

	if (M_export_flag)
	{
		M_Solid_model->computeElementStresses(*solution_rb);
		exporter.postProcess( istance );
		exporter.closeFile();
	}
	// Close HDF5 Files
	HDF5_importer.Close();

	StopTimer("PostProcessing");
	*M_solution = *solution_rb;
	StopGlobalTimer();

}
//=========================================================================
void StructuralSolverTierIII::SolveMDEIM(const std::vector<double>& parameters, meshSubPtr_Type& meshSub, int configuration_number)
{
	StartGlobalTimer();
	//=========================================================================================
	// Hyper-reduction HDF5 file manager
	StartTimer();

	std::string h5_Jac_prefix     = "/"+M_datafile("jacobian_matrix/prefix", "A")+"/";
	std::string h5_resI_prefix    = "/"+M_datafile("residual_internal_force/prefix", "F")+"/";
	std::string h5_resE_prefix    = "/"+M_datafile("residual_external_force/prefix", "F")+"/";

	DenseHDF5 HDF5dense_importer;

	HDF5dense_importer.Open( M_datafile("rom_exporter/hyper_output_file", "ROM_Hyper")+".h5" );

	int offset;
	HDF5dense_importer.ReadIntValue(h5_Jac_prefix, "size_matrix", offset);

	std::vector<long long int> deim_vectorized_GID;
	HDF5dense_importer.ReadVectorLongInt(h5_Jac_prefix, "indices_vec", deim_vectorized_GID);

	std::vector<int> row_index;
	HDF5dense_importer.ReadVectorInt(h5_Jac_prefix, "row_indices", row_index);

	std::vector<int> col_index;
	HDF5dense_importer.ReadVectorInt(h5_Jac_prefix, "col_indices", col_index);

	std::vector<int> deim_dofs_GID;
	HDF5dense_importer.ReadVectorInt("/ReducedMesh/", "dofs", deim_dofs_GID);

	//std::vector<int> volumes_marked;
	//HDF5dense_importer.ReadVectorInt("/ReducedMesh/", "volumes", volumes_marked);

	//std::vector<int> triangles_marked;
	//HDF5dense_importer.ReadVectorInt("/ReducedMesh/", "triangles", triangles_marked);

	Epetra_SerialDenseMatrix PHI_mat;
	HDF5dense_importer.ReadEpetraSerialDenseMatrix(h5_Jac_prefix, "PHI", PHI_mat);

	Epetra_SerialDenseMatrix PHI_rhsI;
	if (M_StructuralModel!="LinearElasticity")
		HDF5dense_importer.ReadEpetraSerialDenseMatrix(h5_resI_prefix, "PHI", PHI_rhsI);

	Epetra_SerialDenseMatrix PHI_rhsE;
	HDF5dense_importer.ReadEpetraSerialDenseMatrix(h5_resE_prefix, "PHI", PHI_rhsE);

	std::vector<int> deim_rhsI_GID;
	if (M_StructuralModel!="LinearElasticity")
		HDF5dense_importer.ReadVectorInt(h5_resI_prefix, "indices", deim_rhsI_GID);

	std::vector<int> deim_rhsE_GID;
	HDF5dense_importer.ReadVectorInt(h5_resE_prefix, "indices", deim_rhsE_GID);

	// Close HDF5 file
	HDF5dense_importer.Close();

	StopTimer("dense HDF5");

	//=========================================================================================
	mapPtr_Type map_rows;
	map_rows.reset ( new MapEpetra ( *M_extractor->getMapUnmarked() ) );

	// VETTORE SOLUZIONE
	vectorPtr_Type solution_rb ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
	solution_rb->zero();

	vectorPtr_Type increment_h ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
	increment_h->zero();

	// Set Exporter
	ExporterHDF5< mesh_Type > exporter(M_datafile, "exporter");
	if (M_export_flag)
	{
		exporter.setMeshProcId( M_mesh, M_comm->MyPID() );
		string FileName =  "SolutionTierIII"+EpetraExt::toString(configuration_number);
		exporter.setPrefix( FileName);
		exporter.setPostDir( "./" );
		exporter.addVariable ( ExporterData< mesh_Type >::VectorField, "SolutionTierIII", M_dispFESpace, solution_rb, UInt ( 0 ) );
	}

	matrixPtr_Type matrix_I_I;
	matrixPtr_Type matrix_I_D;
	vectorPtr_Type rhs_I;
	vectorPtr_Type rhs_D;
	vectorPtr_Type rhs_ext_I;
	vectorPtr_Type rhs_int_I;

	// Pseudo time step for the exporter
	Real istance = 0.0;

	//=============
	// Load Matrix V
	StartTimer();

	// ROM HDF5 file manager
	EpetraExt::HDF5 HDF5_importer(*M_comm);
	HDF5_importer.Open( M_datafile ("importer/file_name", "ROM.h5") );

	int num_basis;
	HDF5_importer.Read("info", "num_basis", num_basis);

	std::vector<int> column_indices;
	for ( int i = 0; i < num_basis; ++i )
		column_indices.push_back(i);

	int* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	boost::shared_ptr<MapEpetra> map_column_V ( new MapEpetra ( -1, static_cast<int> (column_indices.size() ), pointerToDofs, M_comm ) );

	Epetra_Map map_col(*map_column_V->map(Unique));

	Epetra_CrsMatrix* M_V;
	HDF5_importer.Read ( "V", map_col, *M_extractor->getMapUnmarked()->map(Unique), M_V );

	StopTimer("loading matrix V");

	vectorPtr_Type sol ( new vector_Type (*map_column_V, Unique ) );
	vectorPtr_Type increment ( new vector_Type (*map_column_V, Unique ) );
	//=============


	// Jacobian matrix and residual vector
	vectorPtr_Type residual( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_ext( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_int( new vector_Type ( M_dispFESpace->map(), Unique ) );

	// Material model
	if (M_export_flag)
		M_Solid_model->setExportStresses(exporter, M_FESpace_scalar);
	M_Solid_model->updateLoad(M_pressureLoad);

	//=========================================================================================
	// START NEWTON ITERATIONS
	int Newton_MaxIt = M_datafile("newton/maxiter",5);
	int Newton_It = 0;
	double Newton_tolerance = M_datafile("newton/tolerance",1e-6);
	double increment_Norm = Newton_tolerance + 1;
	double residual_Norm = 0;
	double residual_Norm0 = 0;
	double residual_Norm_old = 0;
	bool continue_flag = true;

	M_bc->bcUpdate( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );
	bcManageRhs ( *solution_rb, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );

	//=========================================================================================
	// Assemble External forces once and for all
	M_Solid_model->evaluate_external_forces(solution_rb, residual_ext, meshSub);
	M_extractor->restrictVector( residual_ext, "unmarked", rhs_ext_I );

	DEIM<int> DEIM_rhsE(M_comm, PHI_rhsE, deim_rhsE_GID);
	Epetra_SerialDenseVector rhs_DEIM_E;
	DEIM_rhsE.extractEntriesVector(rhs_ext_I, rhs_DEIM_E);
	Epetra_SerialDenseVector theta_vecE;
	DEIM_rhsE.solveOnline(rhs_DEIM_E, theta_vecE);

	int size_theta_vecE = theta_vecE.Length();
	Epetra_MultiVector resE_n ( *map_column_V->map(Unique), 1 );
	resE_n.PutScalar(0.0);

	Epetra_MultiVector* importedVector;
	for ( int i = 1; i <= size_theta_vecE; ++i )
	{
		HDF5_importer.Read("resEN_"+EpetraExt::toString(i), map_col, importedVector, false );
		resE_n.Update( theta_vecE[i-1], *importedVector, 1.0);
	}

	if ( M_comm->MyPID() == 0 )
		std::cout << "\n\n========= Start Newton Iterations ==========\n\n";

	while (increment_Norm > Newton_tolerance && Newton_It < Newton_MaxIt && continue_flag)
	{
		matrixPtr_Type jacobian ( new matrix_Type ( M_dispETFESpace->map(), 100 ) );

		Newton_It  = Newton_It + 1;
		increment->zero();
		increment_h->zero();
		sol->zero();

		//=========================================================================================
		// ASSEMBLAGGIO SU MESH RIDOTTA
		StartTimer();

		M_Solid_model->evaluate_jacobian(solution_rb, jacobian, meshSub);

		//=========================================================================================
		// ORA DEVO ESTRARRE DALLA MATRICE LE COMPONENTI INTERNE,DI DIRICHLET E MISTE

		M_extractor->restrictMatrix( jacobian, "unmarked", "unmarked", matrix_I_I );
		M_extractor->restrictMatrix( jacobian, "unmarked", "marked",  matrix_I_D );

		// PESCO I VALORI DALLA MATRICE E LI METTO IN VETTORE DENSO
		// DEIM object
		DEIM<long long int> DEIM_mat(M_comm, PHI_mat, deim_vectorized_GID);
		Epetra_SerialDenseVector rhs_MDEIM;
		DEIM_mat.extractEntriesMatrix(jacobian, row_index, col_index, rhs_MDEIM);

		//=========================================================================================
		// solve Interpolation problem Ax = b

		Epetra_SerialDenseVector theta_mat;
		DEIM_mat.solveOnline(rhs_MDEIM, theta_mat);

		StopTimer("assembling matrix + MDEIM");

		//=========================================================================================
		// ASSEMBLAGGIO RHS SU MESH RIDOTTA

		StartTimer();

		M_Solid_model->evaluate_internal_forces(solution_rb, residual_int, meshSub);
		bcManageRhs ( *residual_int, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );
		*residual_int *= -1.0;

		M_extractor->restrictVector( residual_int, "unmarked", rhs_I );
		M_extractor->restrictVector( residual_int, "marked", rhs_D );

		vectorPtr_Type rhs_ridotto ( new VectorEpetra ( *rhs_I - (*matrix_I_D * (*rhs_D) ) ) );

		vectorPtr_Type increment_h_I ( new VectorEpetra ( rhs_I->map() ) );
		increment_h_I->zero();
		Epetra_MultiVector resI_n ( *map_column_V->map(Unique), 1 );

		if (M_StructuralModel!="LinearElasticity")
		{
			DEIM<int> DEIM_rhsI(M_comm, PHI_rhsI, deim_rhsI_GID);
			Epetra_SerialDenseVector rhs_DEIM_I;
			DEIM_rhsI.extractEntriesVector(rhs_ridotto, rhs_DEIM_I);

			Epetra_SerialDenseVector theta_vecI;
			DEIM_rhsI.solveOnline(rhs_DEIM_I, theta_vecI);

			int size_theta_vecI = theta_vecI.Length();
			resI_n.PutScalar(0.0);

			for ( int i = 1; i <= size_theta_vecI; ++i )
			{
				HDF5_importer.Read("resIN_"+EpetraExt::toString(i), map_col, importedVector, false );
				resI_n.Update( theta_vecI[i-1], *importedVector, 1.0);
			}
		}
		else
		{
			resI_n.PutScalar(0.0);
		}
		StopTimer("assembling rhs + DEIM");


		//=========================================================================================
		// Assemble and solve RB problem
		StartTimer();

		Epetra_CrsMatrix * An = new Epetra_CrsMatrix ( Copy, *map_column_V->map ( Unique ), num_basis, false );

		Epetra_CrsMatrix* importedMatrix;

		Amesos_BaseSolver * Solver;
		Amesos Factory;
		std::string SolverType = M_datafile ("solver/method", "Klu");

		Epetra_LinearProblem Problem;

		// form system
		int size_theta_mat = theta_mat.Length();

		for ( int i = 1; i <= size_theta_mat; ++i )
		{
			HDF5_importer.Read("JacN_"+EpetraExt::toString(i), map_col, map_col, importedMatrix );
			EpetraExt::MatrixMatrix::Add ( *importedMatrix, false, theta_mat[i-1], *An, 1. );
		}
		An->FillComplete();

		Epetra_MultiVector fn ( *map_column_V->map(Unique), 1 );
		fn.Update( 1.0,  resI_n, -1.0, resE_n, 0.0);

		fn.Norm2(&residual_Norm);
		if (Newton_It == 1)
			residual_Norm0 = residual_Norm;

		if (Newton_It>3 && residual_Norm > 100*residual_Norm_old)
		{
			if ( M_comm->MyPID() == 0 )
			{
				std::cout << "\n**** Iteration " << Newton_It << ", Residual Norm = " << residual_Norm/residual_Norm0 << "\n";
				std::cout << "===================================================\n";
			}
			continue_flag = false;
		}
		else
		{

			// Only for testing in combination with full assembly
			/*
        	Epetra_FECrsMatrix *An2;

        	boost::shared_ptr<Epetra_FECrsMatrix> tmp ( new Epetra_FECrsMatrix ( Copy, V->OperatorRangeMap(), num_basis ) );
        	int errCode = EpetraExt::MatrixMatrix::Multiply ( *matrix_I_I->matrixPtr(), false, *V, false, *tmp, false );
        	tmp->GlobalAssemble( V->OperatorDomainMap(), V->OperatorRangeMap() );
        	An2 =  new Epetra_FECrsMatrix ( Copy, V->OperatorDomainMap(), num_basis ) ;
        	errCode = EpetraExt::MatrixMatrix::Multiply ( *V, true, *tmp, false, *An2, false );
        	An2->GlobalAssemble(V->OperatorDomainMap(),V->OperatorDomainMap());

        	Epetra_FECrsMatrix *An_err;

        	An_err =  new Epetra_FECrsMatrix ( Copy, V->OperatorDomainMap(), num_basis ) ;
        	errCode = EpetraExt::MatrixMatrix::Multiply ( *V, true, *tmp, false, *An_err, false );
        	An_err->GlobalAssemble(V->OperatorDomainMap(),V->OperatorDomainMap());

        	EpetraExt::MatrixMatrix::Add( *An,  false, -1.0, *An_err, 1.0);
        	double A_n_error = An_err->NormFrobenius();
        	double A_n_norm = An->NormFrobenius();
        	if (M_verbose)
        		std::cout << "\n======== Relative Error A_N = " << A_n_error/A_n_norm << "\n";
			 */

			/*
             Epetra_Vector fn2 ( *map_column_V->map(Unique) );
             V->Multiply(true,  rhs_ridotto->epetraVector(), fn2);// ATTENZIONE

             Epetra_Vector fn3 ( *map_column_V->map(Unique) );
             V->Multiply(true,  rhs_ext_I->epetraVector(), fn3);// ATTENZIONE

             std::cout << resI_n << std::endl;
             std::cout << "\n===============================\n" << std::endl;
             std::cout << fn2 << std::endl;

             Epetra_MultiVector fn4 ( *map_column_V->map(Unique), 1 );
             fn4.Update( 1.0,  fn2, -1.0, fn3, 0.0);
			 */

			// solve
			Problem.SetOperator(An);
			Problem.SetRHS(&fn);
			Problem.SetLHS(&sol->epetraVector());

			Solver = Factory.Create(SolverType, Problem);
			Solver->SymbolicFactorization();
			Solver->NumericFactorization();
			Solver->Solve();

			StopTimer("assembling RB system and solving");

			M_V->Multiply (false, sol->epetraVector(), increment_h_I->epetraVector());
			M_extractor->sumUnmarkedIntoGlobal( increment_h_I, increment_h );

			*solution_rb += *increment_h;

			matrix_I_I->zero();
			matrix_I_D->zero();
			rhs_I->zero();
			rhs_D->zero();
			rhs_ridotto->zero();

			increment_Norm = (increment_h->norm2()) / (solution_rb->norm2());


			if ( M_comm->MyPID() == 0 )
			{
				std::cout << "\n**** Iteration " << Newton_It << ", Norm Increment = " <<  increment_Norm << ", Residual Norm = " << residual_Norm/residual_Norm0 << "\n";
				std::cout << "===================================================\n";
			}

			residual_Norm_old = residual_Norm;

			if (M_StructuralModel=="LinearElasticity")
				continue_flag = false;

			//if (M_export_flag)
			//    {
			//        M_Solid_model->computeElementStresses(*solution_rb);
			//        exporter.postProcess( istance );
			//        istance+=1;
			//    }
		}

	}// END NEWTON ITERATIONS


	//=========================================================================================
	//=========================================================================================
	// Post-processing
	StartTimer();

	if (M_export_flag)
	{
		M_Solid_model->computeElementStresses(*solution_rb);
		exporter.postProcess( istance );
		exporter.closeFile();
	}
	// Close HDF5 Files
	HDF5_importer.Close();

	StopTimer("PostProcessing");
	*M_solution = *solution_rb;
	StopGlobalTimer();

}
//=========================================================================
void StructuralSolverTierIII::StartTimer()
{
	M_tStart = high_resolution_clock::now();
}
//=========================================================================
void StructuralSolverTierIII::StopTimer(const std::string name)
{
	M_tEnd = high_resolution_clock::now();
	if (M_verbose)
		cout << "\n\nElapsed time for " << name << ": " << \
		duration_cast<milliseconds>( M_tEnd - M_tStart ).count() << " ms. \n";

}
//=========================================================================
void StructuralSolverTierIII::StartGlobalTimer()
{
	M_tStartGlobal = high_resolution_clock::now();
}
//=========================================================================
void StructuralSolverTierIII::StopGlobalTimer()
{
	M_tEndGlobal = high_resolution_clock::now();
	if (M_verbose)
		cout << "\n\n**** Elapsed time for Entire Simulation: " << \
		duration_cast<milliseconds>( M_tEndGlobal - M_tStartGlobal ).count() << " ms. ****\n";
}
//=========================================================================
void StructuralSolverTierIII::getSolution(vectorPtr_Type& solution)
{
	solution.reset ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
	*solution = *M_solution;
}
//=========================================================================
} // end namespace LifeV
