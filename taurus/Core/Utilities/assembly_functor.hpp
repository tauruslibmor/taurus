/*
 * lameI functor.hpp
 *
 *  Created on:
 *      Author: davide forti
 */

#ifndef _ASSEMBLYFUNCTOR_HPP_
#define _ASSEMBLYFUNCTOR_HPP_

#include <lifev/core/LifeV.hpp>

using namespace LifeV;

template < typename Return_Type>
class assembly_functor
{

public:
    typedef Return_Type return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {

        Real x = spaceCoordinates[0];
        Real y = spaceCoordinates[1];
        Real z = spaceCoordinates[2];

        return function( 0, x, y, z, 0 );

    }

    assembly_functor ( return_Type (*newFunction)(const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/) )
        : function(newFunction) {  }
    ~assembly_functor() {}

private:

    return_Type (*function)( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/ );

};

#endif
