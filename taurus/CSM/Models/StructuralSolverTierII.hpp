//@HEADER
/*!
    @file
    @author F. Negri
    @date 27-01-2016
*/

#ifndef STRUCTURALSOLVERTIERII_H_
#define STRUCTURALSOLVERTIERII_H_

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

namespace LifeV
{

class StructuralSolverTierII
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

	// constructor
	StructuralSolverTierII (const GetPot& dataFile, const commPtr_Type& communicator);

	// destructor
    ~StructuralSolverTierII() {};

    void Setup(const meshPtr_Type& mesh, bcPtr_Type bc, bcPtr_Type bcMark);//, exporterPtr_Type& exporter);

    void Setup(const meshPtr_Type& mesh, bcPtr_Type bc, bcPtr_Type bcMark, Epetra_CrsMatrix* sol_basis);


    void SetPressureLoad(const double load);

    void SetExportOffset(const int offsetExportedJac, const int offsetExportedRhsI, const int offsetExportedRhsE);

    void GetExportOffset(int& numExportedJac, int& numExportedRhsI, int& numExportedRhsE);

    void Solve(const std::vector<double>& parameters, bool export_matrix_map = false, int configuration_number = 1);

    void WriteOffset();

    boost::shared_ptr<StructuralModel>& getMaterialModelPtr() { return M_Solid_model; }

private:

    void StartTimer();
    void StopTimer(const std::string name);
    void StartGlobalTimer();
    void StopGlobalTimer();


private:

    commPtr_Type M_comm;
    GetPot M_datafile;
    uSpaceStdPtr_Type M_dispFESpace;
    uSpaceStdPtr_Type M_FESpace_scalar;
    uSpaceETAPtr_Type M_dispETFESpace;
    meshPtr_Type M_mesh;
    bcPtr_Type M_bc;
    bcPtr_Type M_bc_markingDOFs;

    boost::shared_ptr<StructuralModel> M_Solid_model;
    boost::shared_ptr<DOF_Extractor> M_extractor;

    double M_pressureLoad;

    std::string M_FEspaceOrder;

    std::string M_NLNSnapshotsCollection;
    std::string M_JacobianApproximation;
    bool M_SplitInternalForces;

    int M_numExportedJac;
    int M_numExportedRhsI;
    int M_numExportedRhsE;
    int M_offsetExportedJac;
    int M_offsetExportedRhsI;
    int M_offsetExportedRhsE;

    std::chrono::high_resolution_clock::time_point M_tStartGlobal;
    std::chrono::high_resolution_clock::time_point M_tEndGlobal;
    std::chrono::high_resolution_clock::time_point M_tStart;
    std::chrono::high_resolution_clock::time_point M_tEnd;

    std::string M_h5_Jac_prefix;
    std::string M_h5_resI_prefix;
    std::string M_h5_resIvol_prefix;
    std::string M_h5_resIiso_prefix;
    std::string M_h5_resE_prefix;

    bool M_verbose;

    Epetra_CrsMatrix* V;
    bool M_flagV;

};

} // namespace LifeV

#endif //  STRUCTURALSOLVERTIERII
