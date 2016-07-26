#include <taurus/CSM/Models/StructuralSolverTierII.hpp>
#include "EpetraExt_Utils.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"

using namespace std;
using namespace std::chrono;

namespace LifeV
{

//=========================================================================
// Constructor //
StructuralSolverTierII::StructuralSolverTierII (const GetPot& dataFile, const commPtr_Type& communicator) :
	M_comm             ( communicator ),
	M_datafile         ( dataFile ),
	M_verbose		   ( false),
	M_numExportedJac   ( 0 ),
	M_numExportedRhsI   ( 0 ),
	M_numExportedRhsE   ( 0 ),
	M_offsetExportedJac   ( 0 ),
	M_offsetExportedRhsI   ( 0 ),
	M_offsetExportedRhsE   ( 0 ),
    M_pressureLoad	      (0.0),
    M_flagV			   ( false )
{
	if ( M_comm->MyPID() == 0 )
	{
		M_verbose = true;
	}
}
//=========================================================================
void StructuralSolverTierII::Setup(const meshPtr_Type& localMeshPtr, bcPtr_Type bc, bcPtr_Type bcMark)//, exporterPtr_Type& exporter)
{

	StartTimer();

	M_bc		     = bc;
	M_bc_markingDOFs = bcMark;
	M_mesh = localMeshPtr;
	//M_exporter = exporter;

	// SPAZI ELEMENTI FINITI
	M_FEspaceOrder = M_datafile ("finite_element/degree", "P1");

	M_dispFESpace.reset (new uSpaceStd_Type (M_mesh, M_FEspaceOrder, 3, M_comm) );
	M_FESpace_scalar.reset (new uSpaceStd_Type (M_mesh, M_FEspaceOrder, 1, M_comm) );
	M_dispETFESpace.reset ( new uSpaceETA_Type ( M_mesh, & ( M_dispFESpace->refFE() ), & ( M_dispFESpace->fe().geoMap() ), M_comm ) );

	// BCHandler to mark dofs
	M_bc_markingDOFs->bcUpdate( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );
	M_bc->bcUpdate( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );

	// Marking dofs
	M_extractor.reset ( new DOF_Extractor ( M_dispFESpace ) ) ;
	M_extractor->setup( M_bc_markingDOFs );

	// hyperreduction settings
	M_NLNSnapshotsCollection   = M_datafile ("hyperreduction/NLN_snapshots_collection", "tierI");
	M_JacobianApproximation    = M_datafile ("hyperreduction/jacobian_approximation", "DEIM");
	M_SplitInternalForces      = M_datafile ("hyperreduction/split_internal_force", false);

	if ( M_verbose )
	{
		std::cout << "\n\n*************** CSM TierII Solver  *********************";
		std::cout << "\n** NonLinear Snapshots Collection   = " << M_NLNSnapshotsCollection;
		std::cout << "\n** Jacobian Approximation Strategy  = " << M_JacobianApproximation;
		std::cout << "\n** Split Internal Forces            = " << std::boolalpha << M_SplitInternalForces;
		std::cout << "\n******************************************************\n\n";
	}

	// Solid Model
	M_Solid_model.reset ( StructuralModel::StructuralModelFactory::instance().createObject (M_datafile("structure/model","none")));
	M_Solid_model->setFESpaces ( M_dispFESpace, M_dispETFESpace );
	M_Solid_model->setUp(M_datafile);

	M_pressureLoad = M_datafile ("structure/pressure_load", 0.0);

	// Set prefix
	M_h5_Jac_prefix     = M_datafile("jacobian_matrix/prefix", "A");
	M_h5_resI_prefix    = M_datafile("residual_internal_force/prefix", "F");
	M_h5_resIiso_prefix = M_datafile("residual_internal_force_iso/prefix", "F_iso");
	M_h5_resIvol_prefix = M_datafile("residual_internal_force_vol/prefix", "F_vol");
	M_h5_resE_prefix    = M_datafile("residual_external_force/prefix", "F_E");
	StopTimer("setup CSM TierII solver");
}
//=========================================================================
void StructuralSolverTierII::Setup(const meshPtr_Type& localMeshPtr, bcPtr_Type bc, bcPtr_Type bcMark, Epetra_CrsMatrix* sol_basis)
{
	Setup(localMeshPtr, bc, bcMark);
	V = sol_basis;
	M_flagV = true;
}
//=========================================================================
void StructuralSolverTierII::SetPressureLoad(const double load)
{
	M_pressureLoad = load;
}
//=========================================================================
void StructuralSolverTierII::Solve(const std::vector<double>& parameters, bool export_matrix_map, int configuration_number )
{
	StartGlobalTimer();

	// Vector used for the simulations
	vectorPtr_Type solution( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type increment( new vector_Type ( M_dispFESpace->map(), Unique ) );
	solution->zero();
	increment->zero();

	// Set Exporter
	ExporterHDF5< mesh_Type > exporter(M_datafile, "exporter");
	exporter.setMeshProcId( M_mesh, M_comm->MyPID() );
	string FileName =  M_datafile ("exporter/name_output_file_tierII", "SolutionTierII")+EpetraExt::toString(configuration_number);
	exporter.setPrefix( FileName);
	exporter.setPostDir( "./" );
	exporter.addVariable ( ExporterData< mesh_Type >::VectorField, "SolutionTierII", M_dispFESpace, solution, UInt ( 0 ) );

	matrixPtr_Type matrix_I_I;
	vectorPtr_Type rhs_I;
	vectorPtr_Type rhs_vol_I;
	vectorPtr_Type rhs_iso_I;
	vectorPtr_Type rhs_D;
	vectorPtr_Type rhs_ext_I;
	vectorPtr_Type rhs_int_I;

	// Pseudo time step for the exporter
	Real istance = 0.0;

	// HDF5 files manager
	EpetraExt::HDF5 HDF5_exporter( *M_comm );
	if (M_NLNSnapshotsCollection=="tierII")
		HDF5_exporter.Open(M_datafile("exporter/name_output_file_system", "SystemSnapshots")+".h5");// assume it already exists


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

	//=============
	// Load Matrix V
	if (!M_flagV)
	{
		StartTimer();
		HDF5_importer.Read ( "V", map_col, *M_extractor->getMapUnmarked()->map(Unique), V );
		StopTimer("loading matrix V");
	}
	HDF5_importer.Close();
	//=============

	// Jacobian matrix and residual vector
	vectorPtr_Type residual( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_ext( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_int( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_int_vol( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_int_iso( new vector_Type ( M_dispFESpace->map(), Unique ) );

	// Material model exporter
	M_Solid_model->setExportStresses(exporter, M_FESpace_scalar);
	M_Solid_model->updateLoad(M_pressureLoad);

	M_bc->bcUpdate( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );

	// Newton iterations
	int Newton_MaxIt = M_datafile("newton/maxiter",5);
	int Newton_It = 0;
	double Newton_tolerance = M_datafile("newton/tolerance",1e-6);
	double increment_Norm = Newton_tolerance + 1;
	double residual_Norm = 0;

	bcManageRhs ( *solution, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );

	StartTimer();
	// Assemble External forces once and for all
	M_Solid_model->evaluate_external_forces(solution, residual_ext);
	//bcManageRhs ( *residual_ext, *displ_FESpace->mesh(), displ_FESpace->dof(), *problemBC, displ_FESpace->feBd(), 1.0, 0.0 );

	M_extractor->restrictVector( residual_ext, "unmarked", rhs_ext_I );
	if (M_NLNSnapshotsCollection=="tierII")
	{
		vectorPtr_Type rhs_ext_full_I ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
		rhs_ext_full_I->zero();
		M_extractor->sumUnmarkedIntoGlobal( rhs_ext_I, rhs_ext_full_I );
		HDF5_exporter.Write(M_h5_resE_prefix+"_"+EpetraExt::toString(M_numExportedRhsE+M_offsetExportedRhsE+1), rhs_ext_full_I->epetraVector());
		M_numExportedRhsE = M_numExportedRhsE + 1;
	}

	residual->zero();
	// Assemble Internal forces
	if (M_SplitInternalForces)
	{
		M_Solid_model->evaluate_internal_forces_iso(solution, residual_int_iso);
		bcManageRhs ( *residual_int_iso, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );//check
		M_Solid_model->evaluate_internal_forces_vol(solution, residual_int_vol);
		bcManageRhs ( *residual_int_vol, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );//check
		*residual += ( (*residual_int_iso) + (*residual_int_vol) + (*residual_ext) );
	}
	else
	{
		M_Solid_model->evaluate_internal_forces(solution, residual_int);
		bcManageRhs ( *residual_int, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );
		*residual += ( (*residual_int) + (*residual_ext) );
	}

	bcManageRhs ( *residual, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );
	double residual_Norm0 = residual->norm2();
	double residual_Norm_old = 1;
	double residual_norm_n_0;
	StopTimer("assembling forces before Newton iterations");

	if ( M_verbose )
			std::cout << "\n\n========= Start Newton Iterations ==========\n\n";

	while (increment_Norm > Newton_tolerance && Newton_It < Newton_MaxIt)
	{
		Newton_It  = Newton_It + 1;

		increment->zero();
		matrixPtr_Type jacobian ( new matrix_Type ( M_dispETFESpace->map(), 100 ) );

		// Jacobian assembly
		StartTimer();
		M_Solid_model->evaluate_jacobian(solution, jacobian);
		StopTimer("jacobian assembly");

		// restrict Jacobian matrix
		M_extractor->restrictMatrix( jacobian, "unmarked", "unmarked", matrix_I_I );

		// write matrix_I_I to h5 file
		if (M_NLNSnapshotsCollection=="tierII" && M_JacobianApproximation=="MDEIM")
		{
			StartTimer();
			HDF5_exporter.Write(M_h5_Jac_prefix+"_"+EpetraExt::toString(M_numExportedJac+M_offsetExportedJac+1), *matrix_I_I->matrixPtr());
			M_numExportedJac = M_numExportedJac + 1;
			StopTimer("export Jacobian matrix");
		}

		vectorPtr_Type rhs_ridotto;
		vectorPtr_Type increment_I;
		if (M_SplitInternalForces)
		{
			*residual_int_iso *= -1.0;
			*residual_int_vol *= -1.0;

			M_extractor->restrictVector( residual_int_iso, "unmarked", rhs_iso_I );
			M_extractor->restrictVector( residual_int_vol, "unmarked", rhs_vol_I );

			rhs_ridotto.reset ( new VectorEpetra ( *rhs_iso_I, Unique ) );

			*rhs_ridotto = (*residual_int_iso) + (*residual_int_vol);
			increment_I.reset ( new VectorEpetra ( rhs_iso_I->map() ) );
		}
		else
		{
			*residual_int *= -1.0;
			M_extractor->restrictVector( residual_int, "unmarked", rhs_I );

			rhs_ridotto.reset ( new VectorEpetra ( *rhs_I, Unique ) );
			increment_I.reset ( new VectorEpetra ( rhs_I->map() ) );
		}
		increment_I->zero();


		// write rhs_ridotto to h5 file
		// vector with noncontiguous map are not correctly written to HDF5
		vectorPtr_Type rhs_full_I ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
		if (M_NLNSnapshotsCollection=="tierII")
		{
			StartTimer();
			if (M_SplitInternalForces)
			{
				rhs_full_I->zero();
				M_extractor->sumUnmarkedIntoGlobal( rhs_iso_I, rhs_full_I );
				HDF5_exporter.Write(M_h5_resIiso_prefix+"_"+EpetraExt::toString(M_numExportedRhsI+M_offsetExportedRhsI+1), rhs_full_I->epetraVector());

				rhs_full_I->zero();
				M_extractor->sumUnmarkedIntoGlobal( rhs_vol_I, rhs_full_I );
				HDF5_exporter.Write(M_h5_resIvol_prefix+"_"+EpetraExt::toString(M_numExportedRhsI+M_offsetExportedRhsI+1), rhs_full_I->epetraVector());
			}
			else
			{
				rhs_full_I->zero();
				M_extractor->sumUnmarkedIntoGlobal( rhs_ridotto, rhs_full_I );
				HDF5_exporter.Write(M_h5_resI_prefix+"_"+EpetraExt::toString(M_numExportedRhsI+M_offsetExportedRhsI+1), rhs_full_I->epetraVector());
			}

			M_numExportedRhsI = M_numExportedRhsI + 1;
			StopTimer("export residual vector");
		}

		*rhs_ridotto -= (*rhs_ext_I);

		///////////////////////////////////////////////
		// projection Matrix
		Epetra_FECrsMatrix *An;

		boost::shared_ptr<Epetra_FECrsMatrix> tmp ( new Epetra_FECrsMatrix ( Copy, V->OperatorRangeMap(), num_basis ) );
		int errCode = EpetraExt::MatrixMatrix::Multiply ( *matrix_I_I->matrixPtr(), false, *V, false, *tmp, false );
		tmp->GlobalAssemble( V->OperatorDomainMap(), V->OperatorRangeMap() );
		An =  new Epetra_FECrsMatrix ( Copy, V->OperatorDomainMap(), num_basis ) ;
		errCode = EpetraExt::MatrixMatrix::Multiply ( *V, true, *tmp, false, *An, false );
		An->GlobalAssemble(V->OperatorDomainMap(),V->OperatorDomainMap());

		// projection Vector
		Epetra_Vector rhs_ridotto_n ( *map_column_V->map(Unique) );
		V->Multiply(true,  rhs_ridotto->epetraVector(), rhs_ridotto_n);

		///////////////////////////////////////////////
		vectorPtr_Type sol ( new vector_Type (*map_column_V, Unique ) );

		Amesos_BaseSolver * Solver;
		Amesos Factory;
		std::string SolverType = "Klu";

		Epetra_LinearProblem Problem;

		Problem.SetOperator(An);
		Problem.SetRHS(&rhs_ridotto_n);
		Problem.SetLHS(&sol->epetraVector());

		Solver = Factory.Create(SolverType, Problem);
		Solver->SymbolicFactorization();
		Solver->NumericFactorization();
		Solver->Solve();

		V->Multiply (false, sol->epetraVector(), increment_I->epetraVector());
		M_extractor->sumUnmarkedIntoGlobal( increment_I, increment );

	    if (Newton_It==1)
	    	rhs_ridotto_n.Norm2(&residual_norm_n_0);


		// -------------------------------------------------------------------
		// start backtracking
		double alpha = 1;
		double t = 1e-4;
		double LinSolvTol = 1e-12;
		vectorPtr_Type solution_new( new vector_Type ( M_dispFESpace->map(), Unique ) );
		solution_new->zero();

		*solution_new = *solution + (alpha * (*increment));

		// Assemble Internal forces
		StartTimer();
		if (M_SplitInternalForces)
		{
			M_Solid_model->evaluate_internal_forces_iso(solution_new, residual_int_iso);
			M_Solid_model->evaluate_internal_forces_vol(solution_new, residual_int_vol);
			*residual_int = (*residual_int_iso) + (*residual_int_vol);
		}
		else
		{
			M_Solid_model->evaluate_internal_forces(solution_new, residual_int);
		}
		StopTimer("internal forces assembly");

		bcManageRhs ( *residual_int, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );

		residual->zero();
		*residual += ( (*residual_int) + (*residual_ext) );
		bcManageRhs ( *residual, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );
		residual_Norm = residual->norm2() / (residual_Norm0);

		int backtrack_iter = 0;
		while (residual_Norm > (1-t*(1-LinSolvTol))*residual_Norm_old && backtrack_iter<0 )
		{
			backtrack_iter = backtrack_iter + 1;


			alpha = alpha * 0.5;
			*solution_new = *solution + (alpha * (*increment));

			// Assemble Internal forces
			StartTimer();
			if (M_SplitInternalForces)
			{
				M_Solid_model->evaluate_internal_forces_iso(solution_new, residual_int_iso);
				M_Solid_model->evaluate_internal_forces_vol(solution_new, residual_int_vol);
				*residual_int = (*residual_int_iso) + (*residual_int_vol);
			}
			else
			{
				M_Solid_model->evaluate_internal_forces(solution_new, residual_int);
			}
			StopTimer("   backtracking: internal forces assembly");

			bcManageRhs ( *residual_int, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );

			residual->zero();
			*residual += ( (*residual_int) + (*residual_ext) );
			bcManageRhs ( *residual, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );
			residual_Norm = residual->norm2() / (residual_Norm0);
			if (M_verbose)
				cout << "\n      Backtracking: Relative residual norm = " << residual_Norm << "\n";
		}
		// -------------------------------------------------------------------


		*increment = alpha * (*increment);
		*solution += *increment;


		increment_Norm = (increment->norm2()) / (solution->norm2());
		residual_Norm = residual->norm2() / (residual_Norm0);;
		residual_Norm_old = residual_Norm;

		matrix_I_I->zero();
		if (M_SplitInternalForces)
		{
			rhs_vol_I->zero();
			rhs_iso_I->zero();
		}
		else
		{
			rhs_I->zero();
		}
		rhs_ridotto->zero();


		double residual_norm_n;
		rhs_ridotto_n.Norm2(&residual_norm_n);
		if ( M_verbose )
		{
			std::cout << "\n\n**** Iteration " << Newton_It << ", Norm Increment = " <<  std::scientific << increment_Norm << \
						 ", Relative Residual Norm = " << residual_Norm << \
						 "\n****   Relative Reduced Residual Norm = " << residual_norm_n/residual_norm_n_0 << "\n";
			std::cout << "===================================================\n\n";
		}

	}

	if (export_matrix_map && M_NLNSnapshotsCollection=="tierII")
	{
		HDF5_exporter.Write(M_h5_Jac_prefix+"_map", *matrix_I_I->matrixPtr());
	}


	HDF5_exporter.Close();

	M_offsetExportedJac  += M_numExportedJac;
	M_offsetExportedRhsI += M_numExportedRhsI;
	M_offsetExportedRhsE += M_numExportedRhsE;

	// Post-processing
	M_Solid_model->computeElementStresses(*solution);
	exporter.postProcess( istance );
	exporter.closeFile();

	StopGlobalTimer();

}
//=========================================================================
void StructuralSolverTierII::WriteOffset()
{
	if (M_NLNSnapshotsCollection=="tierII")
	{
		EpetraExt::HDF5 HDF5_exporter( *M_comm );
		HDF5_exporter.Open(M_datafile("exporter/name_output_file_system", "SystemSnapshots")+".h5");// assume it already exists

		HDF5_exporter.Write("info", "num_"+M_h5_Jac_prefix+"Snapshots", M_offsetExportedJac);
		if (M_SplitInternalForces)
		{
			HDF5_exporter.Write("info", "num_"+M_h5_resIiso_prefix+"Snapshots", M_offsetExportedRhsI);
			HDF5_exporter.Write("info", "num_"+M_h5_resIvol_prefix+"Snapshots", M_offsetExportedRhsI);
		}
		else
		{
			HDF5_exporter.Write("info", "num_"+M_h5_resI_prefix+"Snapshots", M_offsetExportedRhsI);
		}
		HDF5_exporter.Write("info", "num_"+M_h5_resE_prefix+"Snapshots", M_offsetExportedRhsE);

		HDF5_exporter.Close();
	}
}
//=========================================================================
void StructuralSolverTierII::StartTimer()
{
	M_tStart = high_resolution_clock::now();
}
//=========================================================================
void StructuralSolverTierII::StopTimer(const std::string name)
{
	M_tEnd = high_resolution_clock::now();
	if (M_verbose)
		cout << "\n\nElapsed time for " << name << ": " << \
		duration_cast<milliseconds>( M_tEnd - M_tStart ).count()/1000.0 << " s. \n";

}
//=========================================================================
void StructuralSolverTierII::StartGlobalTimer()
{
	M_tStartGlobal = high_resolution_clock::now();
}
//=========================================================================
void StructuralSolverTierII::StopGlobalTimer()
{
	M_tEndGlobal = high_resolution_clock::now();
	if (M_verbose)
			cout << "\n\n**** Elapsed time for Entire Simulation: " << \
			duration_cast<milliseconds>( M_tEndGlobal - M_tStartGlobal ).count()/1000.0 << " s. ****\n";
}
//=========================================================================
void StructuralSolverTierII::SetExportOffset(const int offsetExportedJac, const int offsetExportedRhsI, const int offsetExportedRhsE)
{
	M_offsetExportedJac    = offsetExportedJac;
	M_offsetExportedRhsI   = offsetExportedRhsI;
	M_offsetExportedRhsE   = offsetExportedRhsE;
}
//=========================================================================
void StructuralSolverTierII::GetExportOffset(int& offsetExportedJac, int& offsetExportedRhsI, int& offsetExportedRhsE)
{
	offsetExportedJac  = M_offsetExportedJac;
	offsetExportedRhsI = M_offsetExportedRhsI;
	offsetExportedRhsE = M_offsetExportedRhsE;

}
//=========================================================================
} // end namespace LifeV
