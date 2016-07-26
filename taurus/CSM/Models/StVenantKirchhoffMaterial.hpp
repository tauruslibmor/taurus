#ifndef STVENANTKIRCHHOFFMATERIAL_H
#define STVENANTKIRCHHOFFMATERIAL_H

#include <taurus/CSM/Models/StructuralModel.hpp>

namespace LifeV
{

class StVenantKirchhoffMaterial : public StructuralModel
{
public:

	typedef RegionMesh< LinearTetra > mesh_Type;
    typedef boost::shared_ptr<mesh_Type>  meshPtr_Type;
    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    typedef ETFESpace<mesh_Type, MapEpetra, 3, 3 > ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace_Type> ETFESpacePtr_Type;
    typedef FESpace< mesh_Type, MapEpetra > FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type> FESpacePtr_Type;

    typedef ExpressionDefinitions::deformationGradient_Type tensorF_Type;
    typedef ExpressionDefinitions::determinantTensorF_Type determinantF_Type;
    typedef ExpressionDefinitions::rightCauchyGreenTensor_Type tensorC_Type;
    typedef ExpressionDefinitions::minusTransposedTensor_Type minusT_Type;
    typedef ExpressionDefinitions::traceTensor_Type traceTensor_Type;

    typedef ReducedMesh<mesh_Type> meshSub_Type;
    typedef boost::shared_ptr<meshSub_Type> meshSubPtr_Type;

    StVenantKirchhoffMaterial();

	virtual ~StVenantKirchhoffMaterial();

	void setUp(const GetPot& datafile);

	void updateLoad ( const Real& load) { M_load = load; }

	// Evaluate residual on full mesh
	void evaluate_residual(const vectorPtr_Type& solution, vectorPtr_Type& residual,
						   bool iso_forces=true, bool vol_forces=true, bool external_forces=true );

	// Evaluate residual on reduced mesh
	void evaluate_residual(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub,
						   bool iso_forces=true, bool vol_forces=true, bool external_forces=true);


	void evaluate_internal_forces_iso(const vectorPtr_Type& solution, vectorPtr_Type& residual  )
	{
			evaluate_residual(solution, residual, true, false, false );
	}

	void evaluate_internal_forces_vol(const vectorPtr_Type& solution, vectorPtr_Type& residual  )
	{
			evaluate_residual(solution, residual, false, true, false );
	}

	void evaluate_external_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual )
	{
		evaluate_residual(solution, residual, false, false, true );
	}

	void evaluate_internal_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual  )
	{
		evaluate_residual(solution, residual, true, true, false );
	}


	void evaluate_internal_forces_iso(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub  )
	{
		evaluate_residual(solution, residual, meshSub, true, false, false );
	}

	void evaluate_internal_forces_vol(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub  )
	{
		evaluate_residual(solution, residual, meshSub, false, true, false );
	}

	void evaluate_external_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub )
	{
		evaluate_residual(solution, residual, meshSub, false, false, true );
	}

	void evaluate_internal_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub  )
	{
		evaluate_residual(solution, residual, meshSub, true, true, false );
	}

	// Evaluate jacobian on full mesh
	void evaluate_jacobian(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, bool iso_forces=true, bool vol_forces=true);

	void evaluate_jacobian_iso(const vectorPtr_Type& solution, matrixPtr_Type& jacobian)
	{
		evaluate_jacobian(solution, jacobian, true, false);
	}

	void evaluate_jacobian_vol(const vectorPtr_Type& solution, matrixPtr_Type& jacobian)
	{
			evaluate_jacobian(solution, jacobian, false, true);
	}

	// Evaluate jacobian on reduced mesh
	void evaluate_jacobian(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub, bool iso_forces=true, bool vol_forces=true);

	void evaluate_jacobian_iso(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub)
	{
		evaluate_jacobian(solution, jacobian, meshSub, true, false);
	}

	void evaluate_jacobian_vol(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub)
	{
		evaluate_jacobian(solution, jacobian, meshSub, false, true);
	}

	// Evaluates the stresses
	void computeElementStresses ( const vector_Type& solution);

	void setUseFunctors ( const bool& useFunctors ) { M_useFunctors = useFunctors; }

	void setMu_functor ( const boost::shared_ptr<assembly_functor<Real>>& mufunctor ) { M_evaluationMu = mufunctor; }

	void setLambda_functor ( const boost::shared_ptr<assembly_functor<Real>>& lambdafunctor ) { M_evaluationLambda = lambdafunctor; }

	void updateParameters ( const Real& young, const Real& poisson, const Real& load );

	// void setExportStresses ( ExporterHDF5<mesh_Type> & exporter, const FESpacePtr_Type& fespace );

private:

	Real M_mu;
	Real M_lambda;
	matrixSmall_Type M_identity;
	UInt M_offset;
    Real M_load;
    UInt M_flagLoad;
    bool M_internalPressureLoad;
    bool M_useGravity;
    VectorSmall<3> M_gravity;
    Real M_density;

    // Directional Load
    Real M_directionalLoad;
    UInt M_flagDirectionalLoad;
    bool M_useDirectionalLoad;
    VectorSmall<3> M_directionalLoadVector;

    bool M_useFunctors;
    boost::shared_ptr<assembly_functor<Real>> M_evaluationMu;
    boost::shared_ptr<assembly_functor<Real>> M_evaluationLambda;
};

//! Factory create function
inline StructuralModel * createStVenantKirchhoffMaterial()
{
    return new StVenantKirchhoffMaterial ();
}
namespace
{
static bool S_registerStVenantKirchhoff = StructuralModel::StructuralModelFactory::instance().registerProduct ( "StVenantKirchhoff", &createStVenantKirchhoffMaterial );
}


} // namespace LifeV

#endif // StVenantKirchhoffMaterial_H
