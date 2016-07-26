//@HEADER
/*!
    @file
    @author D. Forti
    @date 10-11-2015
*/

#ifndef MESHCONNECTIVITIES_H_
#define MESHCONNECTIVITIES_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/CurrentFEManifold.hpp>

namespace LifeV
{

class MeshConnectivities
{
public:

	typedef RegionMesh< LinearTetra > mesh_Type;
	typedef FESpace<mesh_Type, MapEpetra>   fespace_Type;
	typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;

	MeshConnectivities ( const fespacePtr_Type& fe_space );

    ~MeshConnectivities() {};

    void getConnectivity(std::vector<std::vector<int>>& volumes, std::vector<std::vector<int>>& triangles );

private:

    fespacePtr_Type M_FESpace;
};

} // namespace LifeV

#endif //  MESHCONNECTIVITIES_H_
