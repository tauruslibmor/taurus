#include <taurus/CSM/Models/StructuralModel.hpp>

using namespace LifeV;

void StructuralModel::setFESpaces ( const FESpacePtr_Type& fespace, const ETFESpacePtr_Type& et_fespace )
{
	M_dispFESpace = fespace;
	M_dispETFESpace = et_fespace;
}

void StructuralModel::setExportStresses ( ExporterHDF5<mesh_Type> & exporter, const FESpacePtr_Type& dFESpace_scalar )
{
	int numVolumes = dFESpace_scalar->mesh()->numVolumes();
	std::vector<int> id_elem;
	for ( int i = 0; i < numVolumes; ++i )
	{
		id_elem.push_back ( dFESpace_scalar->mesh()->element(i).id() );
	}
	int* pointerToDofs (0);
	pointerToDofs = &id_elem[0];
	boost::shared_ptr<MapEpetra> map ( new MapEpetra ( -1, static_cast<int> (id_elem.size() ), pointerToDofs, dFESpace_scalar->map().commPtr() ) );

	M_sigma_xx.reset ( new vector_Type (*map,  Unique ) );
	M_sigma_xy.reset ( new vector_Type (*map,  Unique ) );
	M_sigma_xz.reset ( new vector_Type (*map,  Unique ) );
	M_sigma_yx.reset ( new vector_Type (*map,  Unique ) );
	M_sigma_yy.reset ( new vector_Type (*map,  Unique ) );
	M_sigma_yz.reset ( new vector_Type (*map,  Unique ) );
	M_sigma_zx.reset ( new vector_Type (*map,  Unique ) );
	M_sigma_zy.reset ( new vector_Type (*map,  Unique ) );
	M_sigma_zz.reset ( new vector_Type (*map,  Unique ) );
	M_sigma_vm.reset ( new vector_Type (*map,  Unique ) );

	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sigmaXX", dFESpace_scalar, M_sigma_xx, UInt (0),
			ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sigmaXY", dFESpace_scalar, M_sigma_xy, UInt (0),
			ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sigmaXZ", dFESpace_scalar, M_sigma_xz, UInt (0),
			ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sigmaYX", dFESpace_scalar, M_sigma_yx, UInt (0),
			ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sigmaYY", dFESpace_scalar, M_sigma_yy, UInt (0),
			ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sigmaYZ", dFESpace_scalar, M_sigma_yz, UInt (0),
			ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sigmaZX", dFESpace_scalar, M_sigma_zx, UInt (0),
			ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sigmaZY", dFESpace_scalar, M_sigma_zy, UInt (0),
			ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sigmaZZ", dFESpace_scalar, M_sigma_zz, UInt (0),
			ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "sigmaVonMises", dFESpace_scalar, M_sigma_vm, UInt (0),
				ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );

	M_stresses.reset ( new std::vector<std::vector<vectorPtr_Type>> ());
	M_stresses->resize(3);
	for ( int i = 0; i < 3; ++i)
	{
		(*M_stresses)[i].resize(3);
	}

	(*M_stresses)[0][0] = M_sigma_xx;
	(*M_stresses)[0][1] = M_sigma_xy;
	(*M_stresses)[0][2] = M_sigma_xz;
	(*M_stresses)[1][0] = M_sigma_yx;
	(*M_stresses)[1][1] = M_sigma_yy;
	(*M_stresses)[1][2] = M_sigma_yz;
	(*M_stresses)[2][0] = M_sigma_zx;
	(*M_stresses)[2][1] = M_sigma_zy;
	(*M_stresses)[2][2] = M_sigma_zz;
}
