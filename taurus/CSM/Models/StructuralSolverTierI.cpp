#include <taurus/CSM/Models/StructuralSolverTierI.hpp>
#include <taurus/Core/Utilities/DenseHDF5.hpp>

using namespace std;
using namespace std::chrono;

namespace LifeV
{

//=========================================================================
// Constructor //
StructuralSolverTierI::StructuralSolverTierI (const GetPot& dataFile, const commPtr_Type& communicator) :
	M_comm             ( communicator ),
	M_datafile         ( dataFile ),
	M_verbose		   ( false),
	M_numExportedJac   ( 0 ),
	M_numExportedRhsI   ( 0 ),
	M_numExportedRhsE   ( 0 ),
	M_numExportedSol    ( 0 ),
	M_offsetExportedJac   ( 0 ),
	M_offsetExportedRhsI   ( 0 ),
	M_offsetExportedRhsE   ( 0 ),
	M_offsetExportedSol   ( 0 ),
    M_export_flag         (true),
    M_pressureLoad	      (0.0)
{
	if ( M_comm->MyPID() == 0 )
	{
		M_verbose = true;
	}
}
//=========================================================================
void StructuralSolverTierI::Setup(const meshPtr_Type& localMeshPtr, bcPtr_Type bc, bcPtr_Type bcMark, bool export_flag)//, exporterPtr_Type& exporter)
{
	StartTimer();
	M_bc		     = bc;
	M_bc_markingDOFs = bcMark;
	M_mesh = localMeshPtr;
    M_export_flag = export_flag;

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
	M_NLNSnapshotsCollection   = M_datafile ("hyperreduction/NLN_snapshots_collection", "tierII");
	M_JacobianApproximation    = M_datafile ("hyperreduction/jacobian_approximation", "DEIM");

	if ( M_verbose )
	{
		std::cout << "\n\n*************** CSM TierI Solver  *********************";
		std::cout << "\n** NonLinear Snapshots Collection   = " << M_NLNSnapshotsCollection;
		std::cout << "\n** Jacobian Approximation Strategy  = " << M_JacobianApproximation;
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
	M_h5_resE_prefix    = M_datafile("residual_external_force/prefix", "F");
	M_h5_sol_prefix     = M_datafile("solution/prefix", "sol");

	StopTimer("setup CSM TierI solver");

}
//=========================================================================
void StructuralSolverTierI::SetPressureLoad(const double& load)
{
	M_pressureLoad = load;
}
//=========================================================================
void StructuralSolverTierI::Solve(const std::vector<double>& parameters, bool export_matrix_map, int configuration_number )
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
	string FileName =  M_datafile ("exporter/name_output_file", "SolutionTierI")+EpetraExt::toString(configuration_number);
	exporter.setPrefix( FileName);
	exporter.setPostDir( "./" );
	exporter.addVariable ( ExporterData< mesh_Type >::VectorField, "displacementFEM", M_dispFESpace, solution, UInt ( 0 ) );

	matrixPtr_Type matrix_I_I;
	vectorPtr_Type rhs_I;
	vectorPtr_Type rhs_ext_I;
	vectorPtr_Type rhs_int_I;

	// Pseudo time step for the exporter
	Real istance = 0.0;
	M_numExportedSol = 0;

	// HDF5 files manager
	M_HDF5_exporter.reset (new EpetraExt::HDF5 ( *M_comm ) );
    if (M_export_flag && M_NLNSnapshotsCollection=="tierI")
    	M_HDF5_exporter->Open(M_datafile("exporter/name_output_file_system", "SystemSnapshots")+".h5");// assume it already exists
    
    M_HDF5_exporterSol.reset (new EpetraExt::HDF5 ( *M_comm ) );
    if (M_export_flag)
        M_HDF5_exporterSol->Open(M_datafile("exporter/name_output_file_solution", "Snapshots")+".h5");// assume it already exists

	// Jacobian matrix and residual vector
	vectorPtr_Type residual( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_ext( new vector_Type ( M_dispFESpace->map(), Unique ) );
	vectorPtr_Type residual_int( new vector_Type ( M_dispFESpace->map(), Unique ) );

	// Material model exporter
    M_Solid_model->setExportStresses(exporter, M_FESpace_scalar);

    // set pressure load
	M_Solid_model->updateLoad(M_pressureLoad);

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

	M_extractor->restrictVector( residual_ext, "unmarked", rhs_ext_I );
    if (M_export_flag && M_NLNSnapshotsCollection=="tierI")
    {
        vectorPtr_Type rhs_ext_full_I ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
        rhs_ext_full_I->zero();
        M_extractor->sumUnmarkedIntoGlobal( rhs_ext_I, rhs_ext_full_I );
        M_HDF5_exporter->Write(M_h5_resE_prefix+"_"+EpetraExt::toString(M_numExportedRhsE+M_offsetExportedRhsE+1), rhs_ext_full_I->epetraVector());
        M_numExportedRhsE = M_numExportedRhsE + 1;
    }


	// Assemble Internal forces
	M_Solid_model->evaluate_internal_forces(solution, residual_int);
	bcManageRhs ( *residual_int, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );

	residual->zero();
	*residual += ( (*residual_int) + (*residual_ext) );
	bcManageRhs ( *residual, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );
	double residual_Norm0 = residual->norm2();
	bool residual_norm_flag = true;
	double residual_Norm_old = 1;
	StopTimer("assembling forces before Newton iterations");

	if ( M_verbose )
			std::cout << "\n\n========= Start Newton Iterations ==========\n\n";

	while (increment_Norm > Newton_tolerance && Newton_It < Newton_MaxIt && residual_norm_flag)
	{
		Newton_It  = Newton_It + 1;

		increment->zero();

		// Jacobian assembly
		StartTimer();
		matrixPtr_Type jacobian ( new matrix_Type ( M_dispETFESpace->map(), 100 ) );
		M_Solid_model->evaluate_jacobian(solution, jacobian);
		StopTimer("jacobian assembly");

		// restrict Jacobian matrix
		M_extractor->restrictMatrix( jacobian, "unmarked", "unmarked", matrix_I_I );

		// write matrix_I_I to h5 file
        if (M_export_flag && M_NLNSnapshotsCollection=="tierI" && M_JacobianApproximation=="MDEIM")
        {
            StartTimer();
            M_HDF5_exporter->Write(M_h5_Jac_prefix+"_"+EpetraExt::toString(M_numExportedJac+M_offsetExportedJac+1), *matrix_I_I->matrixPtr());
            M_numExportedJac = M_numExportedJac + 1;
            StopTimer("export Jacobian matrix");
        }

		*residual_int *= -1.0;

		M_extractor->restrictVector( residual_int, "unmarked", rhs_I );

		vectorPtr_Type rhs_ridotto ( new VectorEpetra ( *rhs_I, Unique ) );
		vectorPtr_Type increment_I ( new VectorEpetra ( rhs_I->map() ) );
		increment_I->zero();

		// write rhs_ridotto to h5 file
        if (M_export_flag && M_NLNSnapshotsCollection=="tierI")
        {
            StartTimer();
            vectorPtr_Type rhs_full_I ( new VectorEpetra ( M_dispFESpace->map(), Unique ) );
            rhs_full_I->zero();
            M_extractor->sumUnmarkedIntoGlobal( rhs_ridotto, rhs_full_I );
            M_HDF5_exporter->Write(M_h5_resI_prefix+"_"+EpetraExt::toString(M_numExportedRhsI+M_offsetExportedRhsI+1), rhs_full_I->epetraVector());
            M_numExportedRhsI = M_numExportedRhsI + 1;
            StopTimer("export residual vector");
        }

		*rhs_ridotto -= (*rhs_ext_I);

		// Solve Linear System for the increment
		LinearSolver linearSolver ( M_comm );
		linearSolver.setOperator ( matrix_I_I );

		Teuchos::RCP< Teuchos::ParameterList > aztecList = Teuchos::rcp ( new Teuchos::ParameterList );
		aztecList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

		double LinSolvTol;
		Teuchos::RCP< Teuchos::ParameterList > sublist = Teuchos::rcp ( new Teuchos::ParameterList (aztecList->sublist("Solver: Operator List")));
		Teuchos::RCP< Teuchos::ParameterList > sublist2 = Teuchos::rcp ( new Teuchos::ParameterList (sublist->sublist("Trilinos: AztecOO List")));
		LinSolvTol  = sublist2->get ( "tol", 1e-12 );

		linearSolver.setParameters ( *aztecList );

		prec_Type* precRawPtr;
		basePrecPtr_Type precPtr;
		precRawPtr = new prec_Type;
		precRawPtr->setDataFromGetPot ( M_datafile, "prec" );
		precPtr.reset ( precRawPtr );

		linearSolver.setPreconditioner ( precPtr );

		linearSolver.setRightHandSide( rhs_ridotto );
		linearSolver.solve( increment_I );

		M_extractor->sumUnmarkedIntoGlobal( increment_I, increment );

		// -------------------------------------------------------------------
		// start backtracking
		double alpha = 1;
		double t = 1e-4;
		vectorPtr_Type solution_new( new vector_Type ( M_dispFESpace->map(), Unique ) );
 		solution_new->zero();

		*solution_new = *solution + (alpha * (*increment));

		// Assemble Internal forces
		StartTimer();
		M_Solid_model->evaluate_internal_forces(solution_new, residual_int);
		StopTimer("internal forces assembly");

		bcManageRhs ( *residual_int, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );

		residual->zero();
		*residual += ( (*residual_int) + (*residual_ext) );
		bcManageRhs ( *residual, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_bc, M_dispFESpace->feBd(), 1.0, 0.0 );
		residual_Norm = residual->norm2() / (residual_Norm0);

		int backtrack_iter = 0;
		while (residual_Norm > (1-t*(1-LinSolvTol))*residual_Norm_old && backtrack_iter<1 )
		{
			backtrack_iter = backtrack_iter + 1;

			alpha = alpha * 0.5;
			*solution_new = *solution + (alpha * (*increment));

			// Assemble Internal forces
			StartTimer();
			M_Solid_model->evaluate_internal_forces(solution_new, residual_int);
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
		*solution += (*increment);

		increment_Norm = (increment->norm2()) / (solution->norm2());
		residual_Norm = residual->norm2() / (residual_Norm0);
		residual_Norm_old = residual_Norm;

		if (residual_Norm < 2*LinSolvTol)
			residual_norm_flag = false;

		matrix_I_I->zero();
		rhs_I->zero();
		rhs_ridotto->zero();

		// export to HDF5 the solution
        if (M_export_flag)
        {
            M_HDF5_exporterSol->Write(M_h5_sol_prefix+"_"+EpetraExt::toString(M_numExportedSol+M_offsetExportedSol+1), solution->epetraVector());
            M_numExportedSol += 1;
        }

		if ( M_verbose )
		{
			std::cout << "\n\n**** Iteration " << Newton_It << ", Norm Increment = " <<  std::scientific << increment_Norm << ", Relative Residual Norm = " << residual_Norm << "\n";
			std::cout << "===================================================\n\n";
		}
	}

	if (export_matrix_map && M_export_flag && M_NLNSnapshotsCollection=="tierI" && M_JacobianApproximation=="MDEIM")
	{
		M_HDF5_exporter->Write(M_h5_Jac_prefix+"_map", *matrix_I_I->matrixPtr());
	}

    if (M_export_flag)
    {
    	if (M_export_flag && M_NLNSnapshotsCollection=="tierI")
    		M_HDF5_exporter->Close();

        M_HDF5_exporterSol->Close();
        
        M_offsetExportedSol  += M_numExportedSol;
        M_offsetExportedJac  += M_numExportedJac;
        M_offsetExportedRhsI += M_numExportedRhsI;
        M_offsetExportedRhsE += M_numExportedRhsE;
    }


    // Post-processing
    M_Solid_model->computeElementStresses(*solution);
    exporter.postProcess( istance );
    exporter.closeFile();

    M_solution = solution;
	StopGlobalTimer();

}
//=========================================================================
void StructuralSolverTierI::WriteOffset()
{
	DenseHDF5 HDF5_exporter;
	if (M_NLNSnapshotsCollection=="tierI")
	{
        if (M_verbose)
        {
            HDF5_exporter.Open(M_datafile("exporter/name_output_file_system", "SystemSnapshots")+".h5");// assume it already exists
            
            bool append_system_snapshots = M_datafile("exporter/append_system_snapshots", false);
            if (append_system_snapshots)
            {
                HDF5_exporter.OverWriteIntValue("/info/", "num_"+M_h5_Jac_prefix+"Snapshots", M_offsetExportedJac);
                HDF5_exporter.OverWriteIntValue("/info/", "num_"+M_h5_resI_prefix+"Snapshots", M_offsetExportedRhsI);
                HDF5_exporter.OverWriteIntValue("/info/", "num_"+M_h5_resE_prefix+"Snapshots", M_offsetExportedRhsE);
            }
            else
            {
                HDF5_exporter.WriteIntValue("/info/", "num_"+M_h5_Jac_prefix+"Snapshots", M_offsetExportedJac);
                HDF5_exporter.WriteIntValue("/info/", "num_"+M_h5_resI_prefix+"Snapshots", M_offsetExportedRhsI);
                HDF5_exporter.WriteIntValue("/info/", "num_"+M_h5_resE_prefix+"Snapshots", M_offsetExportedRhsE);
            }
            
            HDF5_exporter.Close();
        }
	}

    if (M_verbose)
    {
        HDF5_exporter.Open(M_datafile("exporter/name_output_file_solution", "SolSnapshots")+".h5");// assume it already exists
        bool append_snapshots = M_datafile("exporter/append_snapshots", false);
        
        if (append_snapshots)
            HDF5_exporter.OverWriteIntValue("/info/", "num_"+M_h5_sol_prefix+"Snapshots", M_offsetExportedSol);
        else
            HDF5_exporter.WriteIntValue("/info/", "num_"+M_h5_sol_prefix+"Snapshots", M_offsetExportedSol);
        
        HDF5_exporter.Close();
    }

}
//=========================================================================
void StructuralSolverTierI::StartTimer()
{
	M_tStart = high_resolution_clock::now();
}
//=========================================================================
void StructuralSolverTierI::StopTimer(const std::string name)
{
	M_tEnd = high_resolution_clock::now();
	if (M_verbose)
		cout << "\n\nElapsed time for " << name << ": " << \
		duration_cast<milliseconds>( M_tEnd - M_tStart ).count()/1000.0 << " s. \n";

}
//=========================================================================
void StructuralSolverTierI::StartGlobalTimer()
{
	M_tStartGlobal = high_resolution_clock::now();
}
//=========================================================================
void StructuralSolverTierI::StopGlobalTimer()
{
	M_tEndGlobal = high_resolution_clock::now();
	if (M_verbose)
			cout << "\n\n**** Elapsed time for Entire Simulation: " << \
			duration_cast<milliseconds>( M_tEndGlobal - M_tStartGlobal ).count()/1000.0 << " s. ****\n";
}
//=========================================================================
void StructuralSolverTierI::SetExportOffset(const int offsetExportedJac, const int offsetExportedRhsI, const int offsetExportedRhsE, const int offsetExportedSol)
{
	M_offsetExportedJac    = offsetExportedJac;
	M_offsetExportedRhsI   = offsetExportedRhsI;
	M_offsetExportedRhsE   = offsetExportedRhsE;
	M_offsetExportedSol    = offsetExportedSol;
}
//=========================================================================
void StructuralSolverTierI::GetExportOffset(int& offsetExportedJac, int& offsetExportedRhsI, int& offsetExportedRhsE, int& offsetExportedSol)
{
	offsetExportedJac  = M_offsetExportedJac;
	offsetExportedRhsI = M_offsetExportedRhsI;
	offsetExportedRhsE = M_offsetExportedRhsE;
	offsetExportedSol  = M_offsetExportedSol;

}
//=========================================================================
void StructuralSolverTierI::getSolution(vectorPtr_Type& solution)
{
    solution = M_solution;
}
//=========================================================================
} // end namespace LifeV
