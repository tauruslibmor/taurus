//@HEADER
/*!
    @file
    @author D. Forti
    @date 10-11-2015
*/

#ifndef IMPORTSNAPSHOTS_H_
#define IMPORTSNAPSHOTS_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <taurus/Core/Utilities/DOF_Extractor.hpp>

namespace LifeV
{

class ImportSnapshots
{
public:

	typedef RegionMesh< LinearTetra > mesh_Type;
	typedef FESpace<mesh_Type, MapEpetra>   fespace_Type;
	typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
	typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<BCHandler> bcPtr_Type;
	typedef boost::shared_ptr<MatrixEpetra<Real>> matrixPtr_Type;

	ImportSnapshots ( const fespacePtr_Type& fe_space, const int numParameters, const GetPot& dataFile, bcPtr_Type bc );

    ~ImportSnapshots() {};

    void readSnapshots(std::string field_name="soluzione");

    void readSnapshotsVectorial(std::string field_name="soluzione");

    void buildColumnMap();

    void getMatrix ( matrixPtr_Type& matrix );

private:

    int M_numParameters;
    GetPot M_datafile;
    fespacePtr_Type M_FESpace;
    mapPtr_Type M_map_column;
    vectorPtr_Type M_solution;
    matrixPtr_Type M_X;
    matrixPtr_Type M_X_restricted;
    boost::shared_ptr<ExporterHDF5<mesh_Type > > M_importer;
    boost::shared_ptr<BCHandler> M_bc_markingDOFs;
    boost::shared_ptr<DOF_Extractor> M_dofExtractor;
};

} // namespace LifeV

#endif //  IMPORTSNAPSHOTS_H_
