#ifndef BCSHEARCUBE_HPP
#define BCSHEARCUBE_HPP

// LifeV includes
#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/BCHandler.hpp>

const int CLAMP   = 3;

namespace LifeV
{

Real fZero (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
	return 0.0;
}

Real markingFunction (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 1.0;
}

typedef boost::shared_ptr<BCHandler> bcPtr_Type;

bcPtr_Type BCh_problem ()
{
	BCFunctionBase zero_function (fZero);
    bcPtr_Type bc (new BCHandler );
    bc->addBC( "CLAMP",  CLAMP,  Essential, Full, zero_function, 3 );
    return bc;
}

bcPtr_Type BCh_marking ()
{
	BCFunctionBase function_marking (markingFunction);
	bcPtr_Type bc (new BCHandler );
    bc->addBC( "CLAMP",  CLAMP,  Essential, Full, function_marking, 3 );
    return bc;
}

}

#endif
