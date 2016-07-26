//@HEADER
/*!
    @file
    @author F. Negri
    @date 02-02-2016
*/

#ifndef BUILDREDUCEDMESH_H_
#define BUILDREDUCEDMESH_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/CurrentFEManifold.hpp>
#include <taurus/Core/HyperReduction/MeshConnectivities.hpp>
#include <taurus/Core/Utilities/DenseHDF5.hpp>


namespace LifeV
{

class BuildReducedMesh
{

public:

	typedef RegionMesh< LinearTetra > mesh_Type;
	typedef FESpace<mesh_Type, MapEpetra>   fespace_Type;
	typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;
	typedef VectorEpetra vector_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<MeshConnectivities> MeshConnectivitiesPtr_Type;

	BuildReducedMesh ( const fespacePtr_Type& fe_space );

    ~BuildReducedMesh() {};

    void appendSelectedDofs(const std::vector<int>& dofs);

    void DofsToNodes_ADR();

    void DofsToNodes_CSM();

    void DofsToNodes_CFD();

    void build( const fespacePtr_Type& FESpace_reducedMesh, vectorPtr_Type&  reduced_mesh, DenseHDF5& HDF5dense_exporter );

    void printDofs();

    void printNodes();

    void setNodes( std::vector<int> vector_nodes );

    void returnVolumes ( std::vector<int>& vector_volumes );

    void setNavierStokes ( const bool& ns );

private:

    void eliminateDuplicates();

    //void DofsToNodes_CFD();

private:

    //fespacePtr_Type M_FESpace_serial;
    int M_numScalarDof;
    std::vector< std::vector<int> > M_mesh_connectivity_volumes;
    std::vector< std::vector<int> > M_mesh_connectivity_triangles;
    MeshConnectivitiesPtr_Type M_connectivity;
    std::vector<int> M_dofs;
    std::vector<int> M_nodes;
    bool M_navierStokes;
};

} // namespace LifeV

#endif //  BUILDREDUCEDMESH_H_
