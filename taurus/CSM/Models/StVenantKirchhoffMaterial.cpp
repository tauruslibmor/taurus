#include <taurus/CSM/Models/StVenantKirchhoffMaterial.hpp>

using namespace std;

using namespace LifeV;

//=============================================================================================
StVenantKirchhoffMaterial::StVenantKirchhoffMaterial():
		M_mu ( 0.0 ),
		M_lambda ( 0.0 ),
		M_offset ( 0 ),
		M_load ( 0.0 ),
		M_flagLoad ( 0 ),
		M_internalPressureLoad ( false ),
		M_useFunctors ( false ),
		M_useDirectionalLoad ( false )
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
StVenantKirchhoffMaterial::~StVenantKirchhoffMaterial()
{}
//=============================================================================================
void StVenantKirchhoffMaterial::setUp(const GetPot& datafile)
{
    // read material parameters
    Real young   = datafile ("structure/young", 1e6 );
    Real poisson = datafile ("structure/poisson", 0.3 );

    M_mu   = young / ( 2.0 * ( 1.0 + poisson ) );
    M_lambda = ( young * poisson ) / ( ( 1.0 + poisson ) * ( 1.0 - 2.0 * poisson ) );

    // read directional load (if it's the case)
    M_useDirectionalLoad = datafile ("structure/useDirectionalLoad", false );

    if ( M_useDirectionalLoad )
    {
    	M_flagDirectionalLoad = datafile ("structure/flagDirectionalLoad", 1 );
    	M_directionalLoad = datafile ("structure/directionalLoad", 0.0 );

    	for ( int i = 0; i < datafile.vector_variable_size ( ("structure/vectorDirectionalLoad") ); ++i )
    	{
    		M_directionalLoadVector(i) = datafile ("structure/vectorDirectionalLoad", 0.0, i);
    	}
    }

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
    	std::cout << "\n** Model: StVenantKirchhoffMaterial";
    	std::cout << "\n** Number of Dofs        = " << M_dispFESpace->dof().numTotalDof()*3;
    	std::cout << "\n** Lame parameter mu     = " << M_mu;
    	std::cout << "\n** Lame parameter lambda = " << M_lambda;
    	std::cout << "\n** InternalPressureLoad  = " << std::boolalpha << M_internalPressureLoad;
    	std::cout << "\n** Gravity               = " << std::boolalpha << M_useGravity;
    	std::cout << "\n******************************************************\n\n";
    }

}
//=============================================================================================
void StVenantKirchhoffMaterial::updateParameters(const Real& young, const Real& poisson, const Real& load)
{
	M_mu   = young / ( 2.0 * ( 1.0 + poisson ) );
	M_lambda = ( young * poisson ) / ( ( 1.0 + poisson ) * ( 1.0 - 2.0 * poisson ) );
	M_directionalLoad = load;
}
//=============================================================================================
void StVenantKirchhoffMaterial::evaluate_residual(const vectorPtr_Type& solution, vectorPtr_Type& residual, bool iso_forces, bool vol_forces, bool external_forces )
{
	residual.reset (new vector_Type ( M_dispFESpace->map() ) );
	residual->zero();

	using namespace ExpressionAssembly;

	vectorPtr_Type solution_rep ( new vector_Type ( *solution, Repeated ) );

	// Definition of F
	tensorF_Type F = ExpressionDefinitions::deformationGradient( M_dispETFESpace, *solution_rep, M_offset, M_identity );

	// Definition of tensor C
	tensorC_Type C = ExpressionDefinitions::tensorC( transpose(F), F );

	// Definition of tr( C )
	traceTensor_Type I_C = ExpressionDefinitions::traceTensor( C );


	if ( iso_forces )
	{

		//Computation of the isochoric part
		if ( M_useFunctors )
		{
			integrate ( elements ( M_dispETFESpace->mesh() ),
					M_dispFESpace->qr(),
					M_dispETFESpace,
					value ( 1.0 / 2.0 ) * eval(M_evaluationLambda, X) * ( I_C - 3.0 ) * dot ( F, grad (phi_i) )

			) >> residual;
		}
		else
		{
			integrate ( elements ( M_dispETFESpace->mesh() ),
					M_dispFESpace->qr(),
					M_dispETFESpace,
					value ( 1.0 / 2.0 ) * value ( M_lambda ) * ( I_C - 3.0 ) * dot ( F, grad (phi_i) )

			) >> residual;
		}

	}

	if ( vol_forces )
	{
		//Computation of the isochoric part
		if ( M_useFunctors )
		{
			integrate ( elements ( M_dispETFESpace->mesh() ),
					M_dispFESpace->qr(),
					M_dispETFESpace,
					value (-1.0) * eval(M_evaluationMu, X) * dot ( F , grad (phi_i) ) +
					eval(M_evaluationMu, X) * dot ( F * C , grad (phi_i) )
			) >> residual;
		}
		else
		{
			integrate ( elements ( M_dispETFESpace->mesh() ),
					M_dispFESpace->qr(),
					M_dispETFESpace,
					value (-1.0) * value ( M_mu ) * dot ( F , grad (phi_i) ) +
					value ( M_mu ) * dot ( F * C , grad (phi_i) )
			) >> residual;
		}
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

		if ( M_useDirectionalLoad )
		{
			QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
			integrate( boundary(M_dispFESpace->mesh(), M_flagDirectionalLoad ),
					myBDQR,
					M_dispETFESpace,
					value(-1.0*M_directionalLoad)* dot(phi_i, M_directionalLoadVector )
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
// Reduced Mesh version
void StVenantKirchhoffMaterial::evaluate_residual(const vectorPtr_Type& solution, vectorPtr_Type& residual, const meshSubPtr_Type& meshSub, bool iso_forces, bool vol_forces, bool external_forces )
{
	residual.reset (new vector_Type ( M_dispFESpace->map() ) );
	residual->zero();

	using namespace ExpressionAssembly;

	vectorPtr_Type solution_rep ( new vector_Type ( *solution, Repeated ) );

	// Definition of F
	tensorF_Type F = ExpressionDefinitions::deformationGradient( M_dispETFESpace, *solution_rep, M_offset, M_identity );

	// Definition of tensor C
	tensorC_Type C = ExpressionDefinitions::tensorC( transpose(F), F );

	// Definition of tr( C )
	traceTensor_Type I_C = ExpressionDefinitions::traceTensor( C );

	if ( iso_forces )
	{
		//Computation of the  isochoric part

		if ( M_useFunctors )
		{
			integrate ( elements (  M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ),
					M_dispFESpace->qr(),
					M_dispETFESpace,
					value ( 1.0 / 2.0 ) * eval(M_evaluationLambda, X) * ( I_C - 3.0 ) * dot ( F, grad (phi_i) )

			) >> residual;
		}
		else
		{
			integrate ( elements (  M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ),
					M_dispFESpace->qr(),
					M_dispETFESpace,
					value ( 1.0 / 2.0 ) * value ( M_lambda ) * ( I_C - 3.0 ) * dot ( F, grad (phi_i) )

			) >> residual;
		}
	}

	if ( vol_forces )
	{
		if ( M_useFunctors )
		{
			integrate ( elements (  M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ),
					M_dispFESpace->qr(),
					M_dispETFESpace,
					value (-1.0) * eval(M_evaluationMu, X) * dot ( F , grad (phi_i) ) +
					eval(M_evaluationMu, X) * dot ( F * C , grad (phi_i) )
			) >> residual;
		}
		else
		{
			integrate ( elements (  M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ),
					M_dispFESpace->qr(),
					M_dispETFESpace,
					value (-1.0) * value ( M_mu ) * dot ( F , grad (phi_i) ) +
					value ( M_mu ) * dot ( F * C , grad (phi_i) )
			) >> residual;
		}
	}

	if ( external_forces )
	{
		if ( M_internalPressureLoad )
		{
			QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );

			integrate( boundary(M_dispFESpace->mesh(), M_flagLoad ),
					myBDQR,
					M_dispETFESpace,
					value(-1.0 * M_load)* dot(phi_i, Nface)
			) >> residual;
		}

		if ( M_useDirectionalLoad )
		{
			QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
			integrate( boundary(M_dispFESpace->mesh(), M_flagDirectionalLoad ),
					myBDQR,
					M_dispETFESpace,
					value(-1.0*M_directionalLoad)* dot(phi_i, M_directionalLoadVector )
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
void StVenantKirchhoffMaterial::evaluate_jacobian(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, bool iso_forces, bool vol_forces )
{
    using namespace ExpressionAssembly;

    jacobian->zero();

    vectorPtr_Type solution_rep ( new vector_Type ( *solution, Repeated ) );

    // Definition of F
    tensorF_Type F = ExpressionDefinitions::deformationGradient( M_dispETFESpace, *solution_rep, M_offset, M_identity );

    // Definition of tensor C
    tensorC_Type C = ExpressionDefinitions::tensorC( transpose(F), F );

    // Definition of tr( C )
    traceTensor_Type I_C = ExpressionDefinitions::traceTensor( C );

    if (iso_forces)
    {
    	if ( M_useFunctors )
    	{
    		integrate ( elements (  M_dispETFESpace->mesh() ) ,
    				M_dispFESpace->qr(),
    				M_dispETFESpace,
    				M_dispETFESpace,
    				eval(M_evaluationLambda, X) * dot ( F, grad (phi_j) ) * dot ( F, grad (phi_i) ) +
    				value ( 1.0 / 2.0 ) * eval(M_evaluationLambda, X) * ( I_C - value (3.0) ) * dot ( grad (phi_j), grad (phi_i) )

    		) >> jacobian;
    	}
    	else
    	{
			integrate ( elements (  M_dispETFESpace->mesh() ) ,
					M_dispFESpace->qr(),
					M_dispETFESpace,
					M_dispETFESpace,
					value ( M_lambda ) * dot ( F, grad (phi_j) ) * dot ( F, grad (phi_i) ) +
					value ( 1.0 / 2.0 ) * value ( M_lambda ) * ( I_C - value (3.0) ) * dot ( grad (phi_j), grad (phi_i) )
			) >> jacobian;
    	}
    }

    if (vol_forces)
    {
    	if ( M_useFunctors )
    	{
    		integrate ( elements (  M_dispETFESpace->mesh() ) ,
    				M_dispFESpace->qr(),
    				M_dispETFESpace,
    				M_dispETFESpace,
    				value (-1.0) * eval(M_evaluationMu, X) * dot ( grad (phi_j), grad (phi_i) ) +
    				eval(M_evaluationMu, X) * dot ( grad (phi_j) * C , grad (phi_i) ) +
    				eval(M_evaluationMu, X) * dot ( F * transpose (grad (phi_j) ) * F  , grad (phi_i) ) +
    				eval(M_evaluationMu, X) *  dot ( F * transpose (F) * grad (phi_j) , grad (phi_i) )
    		) >> jacobian;
    	}
    	else
    	{
    		integrate ( elements (  M_dispETFESpace->mesh() ) ,
    				M_dispFESpace->qr(),
    				M_dispETFESpace,
    				M_dispETFESpace,
    				value (-1.0) * value ( M_mu ) * dot ( grad (phi_j), grad (phi_i) ) +
    				value ( M_mu ) * dot ( grad (phi_j) * C , grad (phi_i) ) +
    				value ( M_mu ) * dot ( F * transpose (grad (phi_j) ) * F  , grad (phi_i) ) +
    				value ( M_mu ) *  dot ( F * transpose (F) * grad (phi_j) , grad (phi_i) )
    		) >> jacobian;
    	}
    }

    jacobian->globalAssemble();
}
//=============================================================================================
// Reduced Mesh version
void StVenantKirchhoffMaterial::evaluate_jacobian(const vectorPtr_Type& solution, matrixPtr_Type& jacobian, const meshSubPtr_Type& meshSub, bool iso_forces, bool vol_forces )
{
    using namespace ExpressionAssembly;

    jacobian->zero();

    vectorPtr_Type solution_rep ( new vector_Type ( *solution, Repeated ) );

    // Definition of F
    tensorF_Type F = ExpressionDefinitions::deformationGradient( M_dispETFESpace, *solution_rep, M_offset, M_identity );

    // Definition of tensor C
    tensorC_Type C = ExpressionDefinitions::tensorC( transpose(F), F );

    // Definition of tr( C )
    traceTensor_Type I_C = ExpressionDefinitions::traceTensor( C );

    if (iso_forces)
    {
    	if ( M_useFunctors )
    	{
    		integrate ( elements (  M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ) ,
    				M_dispFESpace->qr(),
    				M_dispETFESpace,
    				M_dispETFESpace,
    				eval(M_evaluationLambda, X) * dot ( F, grad (phi_j) ) * dot ( F, grad (phi_i) ) +
    				value ( 1.0 / 2.0 ) * eval(M_evaluationLambda, X) * ( I_C - value (3.0) ) * dot ( grad (phi_j), grad (phi_i) )
    		) >> jacobian;
    	}
    	else
    	{
    		integrate ( elements (  M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ) ,
    				M_dispFESpace->qr(),
    				M_dispETFESpace,
    				M_dispETFESpace,
    				value ( M_lambda ) * dot ( F, grad (phi_j) ) * dot ( F, grad (phi_i) ) +
    				value ( 1.0 / 2.0 ) * value ( M_lambda ) * ( I_C - value (3.0) ) * dot ( grad (phi_j), grad (phi_i) )
    		) >> jacobian;
    	}
    }

    if (vol_forces)
    {
    	if ( M_useFunctors )
    	{
    		integrate ( elements (  M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ) ,
    				M_dispFESpace->qr(),
    				M_dispETFESpace,
    				M_dispETFESpace,
    				value (-1.0) * eval(M_evaluationMu, X) * dot ( grad (phi_j), grad (phi_i) ) +
    				eval(M_evaluationMu, X) * dot ( grad (phi_j) * C , grad (phi_i) ) +
    				eval(M_evaluationMu, X) * dot ( F * transpose (grad (phi_j) ) * F  , grad (phi_i) ) +
    				eval(M_evaluationMu, X) *  dot ( F * transpose (F) * grad (phi_j) , grad (phi_i) )
    		) >> jacobian;
    	}
    	else
    	{
    		integrate ( elements (  M_dispETFESpace->mesh(), meshSub->getFlag(), meshSub->getNumElements(), meshSub->getSubmesh(), true ) ,
    				M_dispFESpace->qr(),
    				M_dispETFESpace,
    				M_dispETFESpace,
    				value (-1.0) * value ( M_mu ) * dot ( grad (phi_j), grad (phi_i) ) +
    				value ( M_mu ) * dot ( grad (phi_j) * C , grad (phi_i) ) +
    				value ( M_mu ) * dot ( F * transpose (grad (phi_j) ) * F  , grad (phi_i) ) +
    				value ( M_mu ) *  dot ( F * transpose (F) * grad (phi_j) , grad (phi_i) )
    		) >> jacobian;
    	}
    }

    jacobian->globalAssemble();
}
//=============================================================================================
void StVenantKirchhoffMaterial::computeElementStresses ( const vector_Type& solution )
{
	 QuadratureRule qr ( quadRuleTetra1pt );

	 {
		 using namespace ExpressionAssembly;

		 vectorPtr_Type solution_rep ( new vector_Type ( solution, Repeated ) );

		 // Definition of F
		 tensorF_Type F = ExpressionDefinitions::deformationGradient( M_dispETFESpace, *solution_rep, M_offset, M_identity );

		 // Definition of J
		 determinantF_Type J = ExpressionDefinitions::determinantF( F );

		 // Definition of tensor C
		 tensorC_Type C = ExpressionDefinitions::tensorC( transpose(F), F );

		 // Definition of F^-T
		 minusT_Type  F_T = ExpressionDefinitions::minusT( F );

		 // Definition of tr( C )
		 traceTensor_Type I_C = ExpressionDefinitions::traceTensor( C );

		 if ( M_useFunctors )
		 {
			 ComputeStresses ( elements (  M_dispETFESpace->mesh() ),
					 qr,
					 M_dispETFESpace,
					 value ( 1.0 / 2.0 ) * eval(M_evaluationLambda, X) * ( I_C - 3.0 ) * F +
					 value (-1.0) * eval(M_evaluationMu, X) * F  +
					 eval(M_evaluationMu, X) * F * C
			 ) >> *M_stresses;
		 }
		 else
		 {
			 ComputeStresses ( elements (  M_dispETFESpace->mesh() ),
					 qr,
					 M_dispETFESpace,
					 value ( 1.0 / 2.0 ) * value ( M_lambda ) * ( I_C - 3.0 ) * F +
					 value (-1.0) * value ( M_mu ) * F +
					 value ( M_mu ) * F * C
			 ) >> *M_stresses;
		 }
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
