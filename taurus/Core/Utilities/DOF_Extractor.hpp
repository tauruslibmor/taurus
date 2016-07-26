//@HEADER
/*!
    @file
    @author D. Forti
    @date 10-11-2015
*/

#ifndef DOF_EXTRACTOR_H_
#define DOF_EXTRACTOR_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>

namespace LifeV
{

class DOF_Extractor
{
public:

	typedef RegionMesh< LinearTetra > mesh_Type;
	typedef FESpace<mesh_Type, MapEpetra>   fespace_Type;
	typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
	typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<BCHandler> bcPtr_Type;
	typedef boost::shared_ptr<MatrixEpetra<Real>> matrixPtr_Type;

	DOF_Extractor ( const fespacePtr_Type& fe_space );

    ~DOF_Extractor() {};

    void setup ( bcPtr_Type& bcHandler);

    void setFieldSize( const Real& field_size ) { M_field = field_size; };

    void buildMap ( const std::vector<int> indices, mapPtr_Type& map);

    void restrictMatrix ( const matrixPtr_Type& matrix, const std::string first, const std::string second, matrixPtr_Type& result );

    void restrictSnapshots ( const matrixPtr_Type& matrix, matrixPtr_Type& result );

    void restrictVector ( const vectorPtr_Type& vector, const std::string part, vectorPtr_Type& result );

    void restrictVectorFast ( const vectorPtr_Type& vector, const std::vector<int>& indexes, vectorPtr_Type& result );

    void sumUnmarkedIntoGlobal( const vectorPtr_Type& sol_ridotta, vectorPtr_Type& solutionLapRed );

    void restrictMatrixFast ( const matrixPtr_Type& matrix, const std::vector<int>& indexes, matrixPtr_Type& result);

    mapPtr_Type getMapUnmarked() { return M_mapUnmarkedDofs_loc; }
    
private:

    void fillMaps ( );

    void constructMap ( const mapPtr_Type& input_map, const std::vector<int> listDofs, vectorPtr_Type& numeration, mapPtr_Type& output_map);

    int M_status;
    fespacePtr_Type M_FESpace;
    Real M_dimension;
    Real M_field;

    vectorPtr_Type M_markedDofs;

    std::vector<int> indexes_marked;
    std::vector<int> indexes_unmarked;

    std::vector<int> fullvector_indexes_marked;
    std::vector<int> fullvector_indexes_unmarked;

    mapPtr_Type M_mapMarkedDofs_loc;    // map with holes
    mapPtr_Type M_mapUnmarkedDofs_loc;  // map with holes
};

} // namespace LifeV

#endif //  DOF_EXTRACTOR_H_
