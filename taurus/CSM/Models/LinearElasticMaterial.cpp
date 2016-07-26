#include <taurus/CSM/Models/LinearElasticMaterial.hpp>

using namespace std;

using namespace LifeV;

//=============================================================================================
LinearElasticMaterial::LinearElasticMaterial():
		M_mu ( 0.0 ),
		M_lambda ( 0.0 ),
		M_offset ( 0 ),
		M_load ( 0.0 ),
		M_flagLoad ( 0 ),
		M_internalPressureLoad ( false ),
		M_useGravity ( false )
{
	M_identity (0, 0) = 1.0;
	M_identity (0, 1) = 0.0;
	M_identity (0, 2) = 0.0;
	M_identity (1, 0) = 0.0;
	M_identity (1, 1) = 1.0;
	M_identity (1, 2) = 0.0;
	M_identity (2, 0) = 0.0;
	M_identity (2, 1) = 0.0;
	M_identity (2, 2) = 1.0;
}
//=============================================================================================
LinearElasticMaterial::~LinearElasticMaterial()
{}
//=============================================================================================
void LinearElasticMaterial::setUp(const GetPot& datafile)
{
    // read material parameters
    Real young   = datafile ("structure/young", 1e6 );
    Real poisson = datafile ("structure/poisson", 0.3 );
    
    M_mu   = young / ( 2.0 * ( 1.0 + poisson ) );
    M_lambda = ( young * poisson ) / ( ( 1.0 + poisson ) * ( 1.0 - 2.0 * poisson ) );

    // read pressure load (if it's the case)
    M_internalPressureLoad = datafile ("structure/internalPressureLoad", false );
    
    if (M_internalPressureLoad)
    {
        M_flagLoad = datafile ("structure/flagload", 1 );
    }
    
    M_useGravity = datafile ("structure/applyGravity", false );

    M_density = datafile ("structure/density", 0.0 );

    if ( M_useGravity )
    {
    	for ( int i = 0; i < datafile.vector_variable_size ( ("structure/gravity") ); ++i )
    	{
    		M_gravity(i) = datafile ("structure/gravity", 0.0, i);
    	}
    }

    if ( M_dispFESpace->map().commPtr()->MyPID() == 0 )
    {
    	std::cout << "\n\n*************** Structural Model *********************";
    	std::cout << "\n** Model: LinearElasticMaterial";
    	std::cout << "\n** Number of Dofs        = " <<  M_dispFESpace->dof().numTotalDof()*3;
    	std::cout << "\n** Lame parameter mu     = " << M_mu;
    	std::cout << "\n** Lame parameter lambda = " << M_lambda;
    	std::cout << "\n******************************************************\n\n";
    }

}
//
//(sigmaXX*Normals_X+sigmaXY*Normals_Y+sigmaXZ*Normals_Z)*iHat +
//		(sigmaYX*Normals_X+sigmaYY*Normals_Y+sigmaYZ*Normals_Z)*jHat +
//		(sigmaZX*Normals_X+sigmaZY*Normals_Y+sigmaZZ*Normals_Z)*kHat

//=============================================================================================
void LinearElasticMaterial::evaluate_residual_LE(const vectorPtr_Type& solution, vectorPtr_Type& residual, bool internal_forces, bool external_forces )
{
	residual.reset (new vector_Type ( M_dispFESpace->map() ) );
	residual->zero();

	using namespace ExpressionAssembly;

	if ( internal_forces )
	{
		vectorPtr_Type solution_rep ( new vector_Type ( *solution, Repeated ) );

		integrate ( elements ( M_dispETFESpace->mesh() ),
				M_dispFESpace->qr(),
				M_dispETFESpace,
				value( M_mu ) * dot ( grad(M_dispETFESpace, *solution_rep) + transpose( grad( M_dispETFESpace, *solution_rep ) ), grad(phi_i) ) +
				value( M_lambda ) * trace( grad(M_dispETFESpace, *solution_rep) ) * div(phi_i)
		) >> residual;

	}

	if ( external_forces )
	{
		if ( M_internalPressureLoad )
		{
			QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );

			integrate( boundary(M_dispFESpace->mesh(), M_flagLoad ),
					myBDQR,
					M_dispETFESpace,
					value(-1.0*M_load)* dot(phi_i, Nface)
			) >> residual;
		}

		if ( M_useGravity )
		{
			integrate ( elements ( M_dispETFESpace->mesh() ),
						M_dispFESpace->qr(),
						M_dispETFESpace,
						value(-1.0) * value ( M_density ) * dot ( M_gravity, phi_i)
			) >> residual;
		}
	}

	residual->globalAssemble();

}
//=============================================================================================
void LinearElasticMaterial::evaluate_residual(const vectorPtr_Type& solution, vectorPtr_Type& residual, bool iso_forces, bool vol_forces, bool external_forces )
{
	std::cout << "\n WARNING: Isocoric and volumetric splitting not yet supported for LE material\n";
	evaluate_residual_LE(solution, residual, true, external_forces );
}
//=============================================================================================
// Reduced Mesh version
void LinearElasticMaterial::evaluate_residual_LE(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub,  bool internal_forces, bool external_forces )
{
	residual.reset (new vector_Type ( M_dispFESpace->map() ) );
	residual->zero();

	using namespace ExpressionAssembly;

	if ( internal_forces )
	{
		vectorPtr_Type solution_rep ( new vector_Type ( *solution, Repeated ) );

		integrate ( elements (  M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ),
				M_dispFESpace->qr(),
				M_dispETFESpace,
				value( M_mu ) * dot ( grad(M_dispETFESpace, *solution_rep) + transpose( grad( M_dispETFESpace, *solution_rep ) ), grad(phi_i) ) +
				value( M_lambda ) * trace( grad(M_dispETFESpace, *solution_rep) ) * div(phi_i)
		) >> residual;
	}

	if ( external_forces )
	{
		if ( M_internalPressureLoad )
		{
			QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );

			integrate( boundary(M_dispFESpace->mesh(), M_flagLoad ),
					myBDQR,
					M_dispETFESpace,
					value(-1.0*M_load)* dot(phi_i, Nface)
			) >> residual;
		}

		if ( M_useGravity )
		{
			integrate ( elements ( M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ),
						M_dispFESpace->qr(),
						M_dispETFESpace,
						value(-1.0) * value ( M_density ) * dot ( M_gravity, phi_i)
				) >> residual;
		}
	}
	residual->globalAssemble();

}
//=============================================================================================
void LinearElasticMaterial::evaluate_residual(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub, bool iso_forces, bool vol_forces, bool external_forces )
{
	std::cout << "\n WARNING: Isocoric and volumetric splitting not yet supported for LE material\n";
	evaluate_residual_LE(solution, residual, meshSub, true, external_forces );
}
//=============================================================================================
void LinearElasticMaterial::evaluate_jacobian_LE(const vectorPtr_Type& solution, matrixPtr_Type& jacobian )
{
    using namespace ExpressionAssembly;
    
    jacobian->zero();
    
    //Assembling Volumetric Part
    integrate ( elements (  M_dispETFESpace->mesh() ) ,
               M_dispFESpace->qr(),
               M_dispETFESpace,
               M_dispETFESpace,
               value( M_mu ) * dot ( grad(phi_j) + transpose( grad(phi_j) ), grad(phi_i) ) +
               value( M_lambda ) * div(phi_j) * div(phi_i)
               ) >> jacobian;
    
    jacobian->globalAssemble();
}
//=============================================================================================
void LinearElasticMaterial::evaluate_jacobian(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, bool iso_forces, bool vol_forces)
{
	std::cout << "\n WARNING: Isocoric and volumetric splitting not yet supported for LE material\n";
	evaluate_jacobian_LE(solution, jacobian );
}
//=============================================================================================
// Reduced Mesh version
void LinearElasticMaterial::evaluate_jacobian_LE(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub )
{
    using namespace ExpressionAssembly;

    jacobian->zero();

    //Assembling Volumetric Part
    integrate ( elements (  M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ) ,
               M_dispFESpace->qr(),
               M_dispETFESpace,
               M_dispETFESpace,
               value( M_mu ) * dot ( grad(phi_j) + transpose( grad(phi_j) ), grad(phi_i)) +
               value( M_lambda ) * div(phi_j) * div(phi_i)
               ) >> jacobian;

    jacobian->globalAssemble();
}
//=============================================================================================
void LinearElasticMaterial::evaluate_jacobian(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub, bool iso_forces, bool vol_forces)
{
	std::cout << "\n WARNING: Isocoric and volumetric splitting not yet supported for LE material\n";
	evaluate_jacobian_LE(solution, jacobian, meshSub );
}
//=============================================================================================
void LinearElasticMaterial::computeElementStresses ( const vector_Type& solution )
{
	 QuadratureRule qr ( quadRuleTetra1pt );
	 {
		 using namespace ExpressionAssembly;

		 vectorPtr_Type solution_rep ( new vector_Type ( solution, Repeated ) );

		 ComputeStresses ( elements (  M_dispETFESpace->mesh() ),
				 qr,
				 M_dispETFESpace,
				 value(2.0*M_mu) * value(0.5) * ( grad(M_dispETFESpace, *solution_rep) + transpose( grad( M_dispETFESpace, *solution_rep ) ) )
				 + value(M_lambda) * trace( grad(M_dispETFESpace, *solution_rep) ) * M_identity
		 ) >> *M_stresses;
	 }

	 for ( int i = 0; i < M_dispFESpace->mesh()->numVolumes(); ++i )
	 {
		 int index = M_dispFESpace->mesh()->element(i).id();
		 (*M_sigma_vm)[index] = std::sqrt( pow((*M_sigma_xx)[index] - (*M_sigma_yy)[index],2)
				 	 	 	 	 	 	  +pow((*M_sigma_yy)[index] - (*M_sigma_zz)[index],2)
				 	 	 	 	 	      +pow((*M_sigma_xx)[index] - (*M_sigma_zz)[index],2)
				 	 	 	 	 	      +6*( pow((*M_sigma_yz)[index],2) + pow((*M_sigma_zx)[index],2) + pow((*M_sigma_xy)[index],2)) );
	 }
}
//=============================================================================================
