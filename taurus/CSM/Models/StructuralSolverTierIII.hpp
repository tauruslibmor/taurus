//@HEADER
/*!
    @file
    @author F. Negri
    @date 27-01-2016
*/

#ifndef STRUCTURALSOLVERTIERIII_H_
#define STRUCTURALSOLVERTIERIII_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <chrono>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include "EpetraExt_Utils.h"
#include <taurus/Core/Utilities/DOF_Extractor.hpp>
#include <taurus/CSM/Models/StructuralModel.hpp>
#include <taurus/CSM/Models/LinearElasticMaterial.hpp>
#include <taurus/CSM/Models/StVenantKirchhoffMaterial.hpp>

#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"

#include <taurus/Core/Utilities/DenseHDF5.hpp>
#include <taurus/Core/HyperReduction/DEIM.hpp>


namespace LifeV
{

class StructuralSolverTierIII
{
public:

	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<MatrixEpetra<Real>> matrixPtr_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;

	typedef RegionMesh< LinearTetra > mesh_Type;
	typedef boost::shared_ptr< mesh_Type > meshPtr_Type;
	typedef FESpace< mesh_Type, MapEpetra > uSpaceStd_Type;
	typedef boost::shared_ptr< uSpaceStd_Type > uSpaceStdPtr_Type;
	typedef ETFESpace< mesh_Type, MapEpetra, 3, 3 > uSpaceETA_Type;
	typedef boost::shared_ptr< uSpaceETA_Type > uSpaceETAPtr_Type;
	typedef FESpace<mesh_Type, MapEpetra>::function_Type function_Type;
	typedef MatrixEpetra< Real > matrix_Type;
	typedef VectorEpetra vector_Type;
	typedef LinearSolver::SolverType solver_Type;
	typedef LifeV::Preconditioner basePrec_Type;
	typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;
	typedef PreconditionerIfpack prec_Type;
	typedef boost::shared_ptr<prec_Type> precPtr_Type;
	typedef boost::shared_ptr<BCHandler> bcPtr_Type;
	typedef boost::shared_ptr<ExporterHDF5< mesh_Type >> exporterPtr_Type;
    typedef ReducedMesh<mesh_Type> meshSub_Type;
    typedef boost::shared_ptr<meshSub_Type> meshSubPtr_Type;

    typedef boost::shared_ptr<Epetra_MultiVector> multivecPtr_Type;
    typedef boost::shared_ptr<Epetra_CrsMatrix> crsPtr_Type;

	// constructor
	StructuralSolverTierIII (const GetPot& dataFile, const commPtr_Type& communicator);

	// destructor
    StructuralSolverTierIII() {};

    void Setup(const meshPtr_Type& mesh, bcPtr_Type bc, bcPtr_Type bcMark, bool export_flag = true);
    
    void Setup(const meshPtr_Type& mesh, bcPtr_Type bc, bcPtr_Type bcMark, Epetra_CrsMatrix* sol_basis, bool export_flag = true);

    void Solve(const std::vector<double>& parameters, meshSubPtr_Type& meshSub, int configuration_number = 1);
    
    void getSolution(vectorPtr_Type& solution);

    void setInitialGuess(const vectorPtr_Type& guess);

    double computeFEMResidualNorm();

    void SetPressureLoad(const double load);

    void SetVerbosity(const int& verbosity_level);

    boost::shared_ptr<StructuralModel>& getMaterialModelPtr() { return M_Solid_model; }

private:

    void StartTimer();
    void StopTimer(const std::string name);
    void StartGlobalTimer();
    void StopGlobalTimer();
    
    double computeFEMResidualNorm(const vectorPtr_Type& sol);

    void SolveDEIM(const std::vector<double>& parameters, meshSubPtr_Type& meshSub, int configuration_number = 1);

    void SolveDEIM_split(const std::vector<double>& parameters, meshSubPtr_Type& meshSub, int configuration_number = 1);

    void SolveMDEIM(const std::vector<double>& parameters, meshSubPtr_Type& meshSub, int configuration_number = 1);
    
    void load_V();
    
    void load_HyperData();
    
    void updateExternalForces();

    void compute_residual(const vectorPtr_Type& x, multivecPtr_Type& residual);

    void compute_jacobian_DEIM(const vectorPtr_Type& x, crsPtr_Type& jacobian);
    
    void jacobian_projection(const matrixPtr_Type& jac, const boost::shared_ptr<Epetra_FECrsMatrix>& mask, const Epetra_CrsMatrix* Left, Epetra_FECrsMatrix *& result);

    void backtracing(const vectorPtr_Type& x_k, const vectorPtr_Type& delta, double norm_old, double& alpha);

    //void compute_jacobian_MDEIM(const vectorPtr_Type& solution, crsPtr_Type& jacobian);

private:

    commPtr_Type M_comm;
    GetPot M_datafile;
    uSpaceStdPtr_Type M_dispFESpace;
    uSpaceStdPtr_Type M_FESpace_scalar;
    uSpaceETAPtr_Type M_dispETFESpace;
    meshPtr_Type M_mesh;
    bcPtr_Type M_bc;
    bcPtr_Type M_bc_markingDOFs;
    std::string M_FEspaceOrder;
    boost::shared_ptr<DOF_Extractor> M_extractor;

    boost::shared_ptr<StructuralModel> M_Solid_model;
    double M_pressureLoad;
    bool M_export_flag;
    std::string M_JacobianApproximation;
    std::string M_StructuralModel;
    bool M_SplitInternalForces;
    bool M_verbose;

    boost::shared_ptr<Timer> M_Timer;
    std::chrono::high_resolution_clock::time_point M_tStartGlobal;
    std::chrono::high_resolution_clock::time_point M_tEndGlobal;
    std::chrono::high_resolution_clock::time_point M_tStart;
    std::chrono::high_resolution_clock::time_point M_tEnd;
    
    vectorPtr_Type M_solution;

    boost::shared_ptr<EpetraExt::HDF5> M_HDF5_importer;
    DenseHDF5 M_HDF5dense_importer;
    
    std::string M_h5_resIiso_prefix;
    std::string M_h5_resIvol_prefix;
    std::string M_h5_resE_prefix;
    
    Epetra_CrsMatrix* M_V;
    bool M_flagV;
    mapPtr_Type M_map_column_V;
    boost::shared_ptr<Epetra_Map> M_map_col;
    int M_num_basis;

    vectorPtr_Type M_NewtonGuess;
    bool M_flagGuess;
    int M_backtracingIter;
    double M_backtracingFactor;
    double M_residual_Norm0;
    
    // DEIM stuff
    std::vector<int> M_deim_dofs_GID;
    std::vector<int> M_deim_dofs_attached;
    Epetra_SerialDenseMatrix M_PHI_rhsIiso;
    Epetra_SerialDenseMatrix M_PHI_rhsIvol;
    Epetra_SerialDenseMatrix M_PHI_rhsE;
    std::vector<int> M_deim_rhsIiso_GID;
    std::vector<int> M_deim_rhsIvol_GID;
    std::vector<int> M_deim_rhsE_GID;
    Epetra_CrsMatrix* M_LV_iso;
    Epetra_CrsMatrix* M_LV_vol;
    
    boost::shared_ptr<Epetra_Map> M_map_col_LV_iso;
    boost::shared_ptr<Epetra_Map> M_map_col_LV_vol;
    boost::shared_ptr<Epetra_FECrsMatrix> P_DEIM_iso;
    boost::shared_ptr<Epetra_FECrsMatrix> P_DEIM_vol;
    std::vector<Epetra_MultiVector> M_resIiso_n_imported;
    std::vector<Epetra_MultiVector> M_resIvol_n_imported;
    std::vector<Epetra_MultiVector> M_resE_n_imported;

    multivecPtr_Type M_resE_n;
    multivecPtr_Type M_resIiso_n;
    multivecPtr_Type M_resIvol_n;
    
    meshSubPtr_Type M_meshSub;
};

} // namespace LifeV

#endif //  STRUCTURALSOLVERTIERIII
