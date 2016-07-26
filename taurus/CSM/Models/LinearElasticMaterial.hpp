#ifndef LINEARELASTICMATERIAL_H
#define LINEARELASTICMATERIAL_H

#include <taurus/CSM/Models/StructuralModel.hpp>

namespace LifeV
{

class LinearElasticMaterial : public StructuralModel
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

    LinearElasticMaterial();

	virtual ~LinearElasticMaterial();

	void setUp(const GetPot& datafile);

	void updateLoad ( const Real& load) { M_load = load; }

	// Evaluate residual on full mesh
	void evaluate_residual(const vectorPtr_Type& solution, vectorPtr_Type& residual, bool iso_forces=true, bool vol_forces=true, bool external_forces=true );

	void evaluate_residual_LE(const vectorPtr_Type& solution, vectorPtr_Type& residual, bool internal_forces=true, bool external_forces=true );

	// Evaluate residual on reduced mesh
	void evaluate_residual(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub, bool iso_forces=true, bool vol_forces=true, bool external_forces=true );

	void evaluate_residual_LE(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub, bool internal_forces=true, bool external_forces=true);


	void evaluate_external_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual )
	{
		evaluate_residual(solution, residual, false, true );
	}

	void evaluate_internal_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual  )
	{
		evaluate_residual(solution, residual, true, false );
	}


	void evaluate_external_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub )
	{
		evaluate_residual(solution, residual, meshSub, false, true );
	}

	void evaluate_internal_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub  )
	{
		evaluate_residual(solution, residual, meshSub, true, false );
	}

	// Evaluate jacobian on full mesh
	void evaluate_jacobian(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, bool iso_forces=true, bool vol_forces=true);

	void evaluate_jacobian_LE(const vectorPtr_Type& solution, matrixPtr_Type& jacobian);

	// Evaluate jacobian on reduced mesh
	void evaluate_jacobian(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub, bool iso_forces=true, bool vol_forces=true);

	void evaluate_jacobian_LE(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub);

	// Evaluates the stresses
	void computeElementStresses ( const vector_Type& solution);

	// void setExportStresses ( ExporterHDF5<mesh_Type> & exporter, const FESpacePtr_Type& fespace );

private:

	Real M_lambda;
	Real M_mu;
	matrixSmall_Type M_identity;
	UInt M_offset;
    Real M_load;
    UInt M_flagLoad;
    bool M_internalPressureLoad;
    bool M_useGravity;
    VectorSmall<3> M_gravity;
    Real M_density;

};

//! Factory create function
inline StructuralModel * createLinearElasticMaterial()
{
    return new LinearElasticMaterial ();
}
namespace
{
static bool S_registerLinearElasticity = StructuralModel::StructuralModelFactory::instance().registerProduct ( "LinearElasticity", &createLinearElasticMaterial );
}


} // namespace LifeV

#endif // LINEARELASTICMATERIAL_H
