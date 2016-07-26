#ifndef STRUCTURALMODEL_HPP
#define STRUCTURALMODEL_HPP

#include "Epetra_MpiComm.h"

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/Assembly.hpp>
#include <taurus/CSM/Models/AssemblyElementalStructure.hpp>
#include <taurus/CSM/Models/ExpressionDefinitions.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <taurus/Core/HyperReduction/ReducedMesh.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/util/Factory.hpp>
#include <taurus/Core/Utilities/assembly_functor.hpp>

namespace LifeV
{

class StructuralModel
{
public:

	typedef RegionMesh< LinearTetra > mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<vector_Type > vectorPtr_Type;
    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    typedef MapEpetra map_Type;
    typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
    typedef FactorySingleton<Factory<StructuralModel, std::string> > StructuralModelFactory;
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

    StructuralModel(){};

    virtual ~StructuralModel() {}

    void setFESpaces ( const FESpacePtr_Type& fespace, const ETFESpacePtr_Type& et_fespace );

    virtual void setUp(const GetPot& datafile) = 0;

    // Evaluate residual on full mesh
    virtual void evaluate_residual(const vectorPtr_Type& solution, vectorPtr_Type& residual, bool iso_forces=true, bool vol_forces=true, bool external_forces=true ) = 0;

    // Evaluate residual on reduced mesh
    virtual void evaluate_residual(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub,bool iso_forces=true, bool vol_forces=true, bool external_forces=true) = 0;

    virtual void evaluate_external_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual ) = 0;

    virtual void evaluate_internal_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual  ) = 0;

    virtual void evaluate_external_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub ) = 0;

    virtual void evaluate_internal_forces(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub  ) = 0;

    virtual void evaluate_internal_forces_iso(const vectorPtr_Type& solution, vectorPtr_Type& residual  ) {};

    virtual void evaluate_internal_forces_vol(const vectorPtr_Type& solution, vectorPtr_Type& residual  ) {};

    virtual void evaluate_internal_forces_iso(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub  ){};

    virtual void evaluate_internal_forces_vol(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub  ){};


    // Evaluate jacobian on full mesh
    virtual void evaluate_jacobian(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, bool iso_forces=true, bool vol_forces=true) = 0;

    virtual void evaluate_jacobian_iso(const vectorPtr_Type& solution, matrixPtr_Type& jacobian){};

    virtual void evaluate_jacobian_vol(const vectorPtr_Type& solution, matrixPtr_Type& jacobian){};

    // Evaluate jacobian on reduced mesh
    virtual void evaluate_jacobian(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub, bool iso_forces=true, bool vol_forces=true) = 0 ;

    virtual void evaluate_jacobian_iso(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub){};

    virtual void evaluate_jacobian_vol(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub){};

    // Evaluates the element stresses
    virtual void computeElementStresses ( const vector_Type& solution) = 0;

    // Evaluates the element stresses
    virtual void updateLoad ( const Real& /*load*/) {};

    // Set exporter for stresses
    void setExportStresses ( ExporterHDF5<mesh_Type> & exporter, const FESpacePtr_Type& fespace );

    virtual void setUseFunctors ( const bool& /*useFunctors*/ ) { };

    virtual void setMu_functor ( const boost::shared_ptr<assembly_functor<Real>>& /*mufunctor*/ ) { };

    virtual void setBulk_functor ( const boost::shared_ptr<assembly_functor<Real>>& /*bulkfunctor*/ ) { };

    virtual void setLambda_functor ( const boost::shared_ptr<assembly_functor<Real>>& /*lambdafunctor*/ ) { };

    virtual void test( vector_Type& /*residual*/ ) { } ;

    virtual void updateParameters ( const Real& /*young*/, const Real& /*poisson*/, const Real& /*load*/ ) { };

protected:

    FESpacePtr_Type M_dispFESpace;
    ETFESpacePtr_Type M_dispETFESpace;

    vectorPtr_Type M_sigma_xx;
    vectorPtr_Type M_sigma_xy;
    vectorPtr_Type M_sigma_xz;
    vectorPtr_Type M_sigma_yx;
    vectorPtr_Type M_sigma_yy;
    vectorPtr_Type M_sigma_yz;
    vectorPtr_Type M_sigma_zx;
    vectorPtr_Type M_sigma_zy;
    vectorPtr_Type M_sigma_zz;
    vectorPtr_Type M_sigma_vm;

    boost::shared_ptr<std::vector<std::vector<vectorPtr_Type>>> M_stresses;
};

}

#endif // STRUCTURALMODEL_HPP
