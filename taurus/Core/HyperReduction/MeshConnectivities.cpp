#include <taurus/Core/HyperReduction/MeshConnectivities.hpp>

namespace LifeV
{


MeshConnectivities::MeshConnectivities ( const fespacePtr_Type& fe_space ) :
	M_FESpace    ( fe_space )
{
}

void
MeshConnectivities::getConnectivity ( std::vector<std::vector<int>>& volumes, std::vector<std::vector<int>>& triangles  )
{
	// --------------------------------------------
	// CONNETTIVITA' SU VOLUMI
	// --------------------------------------------

	int numTotalDof = M_FESpace->dof().numTotalDof();
	int numVolumes = M_FESpace->mesh()->numVolumes();
	UInt numberLocalDof (M_FESpace->dof().numLocalDof() );

	for ( int i = 0; i < numVolumes; ++i )
	{
		for ( int j = 0; j < numberLocalDof; ++j )
		{
			ID globalDofID (M_FESpace->dof().localToGlobalMap (i, j) );
			volumes[globalDofID].push_back(i);
		}
	}

	// --------------------------------------------
	// CONNETTIVITA' SU TRIANGOLI
	// --------------------------------------------

	CurrentFEManifold feBd1 ( M_FESpace->refFE().boundaryFE(), getGeometricMap ( *M_FESpace->mesh() ).boundaryMap() );

	int numBFaces  = M_FESpace->mesh()->numBFaces();

	for ( int i = 0; i < numBFaces; ++i )
	{
		feBd1.update ( M_FESpace->mesh()->boundaryFacet ( i ), UPDATE_ONLY_CELL_NODES );  // Updating facet information on mesh1
		std::vector<ID> localToGlobalMapOnBFacet1 = M_FESpace->dof().localToGlobalMapOnBdFacet (i);
		for (ID lDof1 = 0; lDof1 < localToGlobalMapOnBFacet1.size(); lDof1++)
		{
			triangles[localToGlobalMapOnBFacet1[lDof1]].push_back(i);
		}
	}
}

} // end namespace LifeV
