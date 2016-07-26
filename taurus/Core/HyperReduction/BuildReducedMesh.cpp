#include <taurus/Core/HyperReduction/BuildReducedMesh.hpp>
#include <taurus/Core/Utilities/Timer.hpp>

namespace LifeV
{

//=========================================================================
// Constructor //
BuildReducedMesh::BuildReducedMesh ( const fespacePtr_Type& fe_space ):
		M_navierStokes ( false )
{
	M_numScalarDof = fe_space->dof().numTotalDof();
	M_mesh_connectivity_volumes.resize( M_numScalarDof );
	M_mesh_connectivity_triangles.resize( M_numScalarDof );
	M_connectivity.reset ( new MeshConnectivities ( fe_space ) );// check!!!!
	M_connectivity->getConnectivity(M_mesh_connectivity_volumes, M_mesh_connectivity_triangles);
}
//=========================================================================
void
BuildReducedMesh::setNavierStokes ( const bool& ns )
{
	M_navierStokes = ns;
}
//=========================================================================
void BuildReducedMesh::appendSelectedDofs ( const std::vector<int>& dofs  )
{

	for ( int k = 0; k < dofs.size(); ++k )
	{
		M_dofs.push_back(dofs[k]);
	}

	eliminateDuplicates();

}
//=========================================================================
void BuildReducedMesh::DofsToNodes_ADR (  )
{

	for ( int k = 0; k < M_dofs.size(); ++k )
	{
		M_nodes.push_back(M_dofs[k]);
	}

}
//=========================================================================
void BuildReducedMesh::DofsToNodes_CSM (  )
{

	for ( int k = 0; k < M_dofs.size(); ++k )
	{
		int c = (int)(M_dofs[k]/M_numScalarDof);
		int r = M_dofs[k] - (c)*M_numScalarDof;
		M_nodes.push_back(r);
	}

	sort( M_nodes.begin(), M_nodes.end() );
	M_nodes.erase( unique( M_nodes.begin(), M_nodes.end() ), M_nodes.end() );

}
//=========================================================================
void BuildReducedMesh::DofsToNodes_CFD (  )
{

	for ( int k = 0; k < M_dofs.size(); ++k )
	{
		int c = (int)(M_dofs[k]/M_numScalarDof);
		int r = M_dofs[k] - (c)*M_numScalarDof;
		M_nodes.push_back(r);
	}

	sort( M_nodes.begin(), M_nodes.end() );
	M_nodes.erase( unique( M_nodes.begin(), M_nodes.end() ), M_nodes.end() );

}
//=========================================================================
void BuildReducedMesh::eliminateDuplicates (  )
{

	sort( M_dofs.begin(), M_dofs.end() );
	M_dofs.erase( unique( M_dofs.begin(), M_dofs.end() ), M_dofs.end() );


}
//=========================================================================
void BuildReducedMesh::setNodes ( std::vector<int> vector_nodes )
{
	M_nodes = vector_nodes;
}
//=========================================================================
void BuildReducedMesh::returnVolumes ( std::vector<int>& vector_volumes )
{
	for ( int k = 0; k < M_nodes.size(); ++k )
	{
		for ( int z = 0; z < M_mesh_connectivity_volumes[M_nodes[k]].size(); ++z )
		{
			vector_volumes.push_back(M_mesh_connectivity_volumes[M_nodes[k]][z]);
		}
	}

	sort( vector_volumes.begin(), vector_volumes.end() );
	vector_volumes.erase( unique( vector_volumes.begin(), vector_volumes.end() ), vector_volumes.end() );
}
//=========================================================================
void BuildReducedMesh::build ( const fespacePtr_Type& FESpace_reducedMesh, vectorPtr_Type&  reduced_mesh, DenseHDF5& HDF5dense_exporter)
{

	Timer myTimer(FESpace_reducedMesh->map().commPtr(), 1.0, "s");
	myTimer.StartTimer();

	//	VETTORE CHE CONTIENE ELEMENTI SU CUI SI DOVRA ASSEMBLARE MESH RIDOTTA
	std::vector<int> volumes_marked;
	std::vector<int> triangles_marked;

	for ( int k = 0; k < M_nodes.size(); ++k )
	{
		for ( int z = 0; z < M_mesh_connectivity_volumes[M_nodes[k]].size(); ++z )
		{
			volumes_marked.push_back(M_mesh_connectivity_volumes[M_nodes[k]][z]);
		}

		for ( int z = 0; z < M_mesh_connectivity_triangles[M_nodes[k]].size(); ++z )
		{
			triangles_marked.push_back(M_mesh_connectivity_triangles[M_nodes[k]][z]);
		}
	}

	sort( volumes_marked.begin(), volumes_marked.end() );
	volumes_marked.erase( unique( volumes_marked.begin(), volumes_marked.end() ), volumes_marked.end() );

	sort( triangles_marked.begin(), triangles_marked.end() );
	triangles_marked.erase( unique( triangles_marked.begin(), triangles_marked.end() ), triangles_marked.end() );

//	if ( FESpace_reducedMesh->map().commPtr()->MyPID() == 0 )
//	{
//		std::cout << "\n============================\n";
//		std::cout << "Volumes for reduced assembly GID:\n\n";
//		for ( int k = 0; k < volumes_marked.size(); ++k )
//		{
//			std::cout << volumes_marked[k] << " ";
//		}
//		std::cout << "\n============================\n";
//	}
//
//	if ( FESpace_reducedMesh->map().commPtr()->MyPID() == 0 )
//	{
//		std::cout << "\n============================\n";
//		std::cout << "Triangles for reduced assembly GID:\n\n";
//		for ( int k = 0; k < triangles_marked.size(); ++k )
//		{
//			std::cout << triangles_marked[k] << " ";
//		}
//		std::cout << "\n============================";
//	}

	// SAVE REDUCED MESH STRUCTURES TO HDF5
	if ( FESpace_reducedMesh->map().commPtr()->MyPID() == 0 )
	{
		HDF5dense_exporter.WriteVectorInt("/ReducedMesh/", "dofs", M_dofs);
		HDF5dense_exporter.WriteVectorInt("/ReducedMesh/", "nodes", M_nodes);
		HDF5dense_exporter.WriteVectorInt("/ReducedMesh/", "volumes", volumes_marked);
		HDF5dense_exporter.WriteVectorInt("/ReducedMesh/", "triangles", triangles_marked);
	}

	// LOOP OVER MARKED VOLUMES TO FILL VECTOR TO VISUALIZE REDUCED MESH
	// and GENERATE LIST OF ALL ASSOCIATED DOFS
	std::vector<int> dofs_attached;

	for ( int k = 0; k < volumes_marked.size(); ++k )
	{
		for ( int z = 0; z < FESpace_reducedMesh->dof().numLocalDof(); ++z)
		{
			int globalDofID (FESpace_reducedMesh->dof().localToGlobalMap (volumes_marked[k], z) );
			(*reduced_mesh)[globalDofID] = 1.0;

			for ( int d = 0; d < nDimensions; ++d )
				dofs_attached.push_back(globalDofID+d*M_numScalarDof);
		}
		if ( M_navierStokes )
		{
			for ( int z = 0; z < 4; ++z)
			{
				int globalDofID (FESpace_reducedMesh->dof().localToGlobalMap (volumes_marked[k], z) );
				dofs_attached.push_back(globalDofID+3*M_numScalarDof);
			}
		}
	}

	CurrentFEManifold feBd1 ( FESpace_reducedMesh->refFE().boundaryFE(), getGeometricMap ( *FESpace_reducedMesh->mesh() ).boundaryMap() );

	for ( int k = 0; k < triangles_marked.size(); ++k )
	{
		feBd1.update ( FESpace_reducedMesh->mesh()->boundaryFacet ( triangles_marked[k] ), UPDATE_ONLY_CELL_NODES );
		std::vector<ID> localToGlobalMapOnBFacet1 = FESpace_reducedMesh->dof().localToGlobalMapOnBdFacet (triangles_marked[k]);

		for ( int z = 0; z < localToGlobalMapOnBFacet1.size(); ++z)
		{
			int globalDofID = (int)(localToGlobalMapOnBFacet1[z]);

			for ( int d = 0; d < nDimensions; ++d )
				dofs_attached.push_back(globalDofID+d*M_numScalarDof);

			if ((*reduced_mesh)[globalDofID] == 1.0)
				(*reduced_mesh)[globalDofID] = 2.0;
			else
				(*reduced_mesh)[globalDofID] = 3.0;
		}

		if ( M_navierStokes )
		{
			for ( int z = 0; z < 3; ++z)
			{
				int globalDofID = (int)(localToGlobalMapOnBFacet1[z]);
				dofs_attached.push_back(globalDofID+3*M_numScalarDof);
			}
		}
	}

	sort( dofs_attached.begin(), dofs_attached.end() );
	dofs_attached.erase( unique( dofs_attached.begin(), dofs_attached.end() ), dofs_attached.end() );
	if ( FESpace_reducedMesh->map().commPtr()->MyPID() == 0 )
	{
		HDF5dense_exporter.WriteVectorInt("/ReducedMesh/", "dofs_attached", dofs_attached);
	}

	myTimer.StopTimer("generating and saving the Reduced Mesh");
	if ( FESpace_reducedMesh->map().commPtr()->MyPID() == 0 )
		std::cout << "============================\n";

}
//=========================================================================
} // end namespace LifeV
