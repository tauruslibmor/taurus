//@HEADER
/*!
    @file
    @author D. Forti
    @date 10-11-2015
*/

#ifndef IMPORTMATRIXSNAPSHOTS_H_
#define IMPORTMATRIXSNAPSHOTS_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
//#include <taurus/Core/Utilities/DOF_Extractor.hpp>
#include "EpetraExt_Utils.h"
#include <EpetraExt_Transpose_RowMatrix.h>

namespace LifeV
{

class ImportMatrixSnapshots
{
public:

	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;
	//typedef RegionMesh< LinearTetra > mesh_Type;

	typedef MapEpetra map_Type;
	typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
    typedef boost::shared_ptr<Epetra_Map> mapEpetraPtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<MatrixEpetra<Real>> matrixPtr_Type;

	ImportMatrixSnapshots ( const mapPtr_Type& map_rows, const commPtr_Type& communicator, const int numParameters, const GetPot& dataFile, const std::string& dataFile_section );

    ~ImportMatrixSnapshots() {};

    void readSnapshots();

    void buildColumnMap();

    void getMatrix ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix );

    void buildRowMap ( const boost::shared_ptr<Epetra_CrsMatrix>& matrix );

    //void buildMatrixRowMap ( );

    int getOffset() const { return M_offset; };

    mapEpetraPtr_Type getRowMap_FEmatrix() { return M_map_matrix_rows->map(Unique); };
    
private:

    void transpose( Epetra_CrsMatrix* input_matrix);

    
    int M_numParameters;
    GetPot M_datafile;
    mapEpetraPtr_Type M_map_column;
    mapEpetraPtr_Type M_map_rows;
    mapPtr_Type M_map_matrix_rows;
    commPtr_Type M_comm;
    std::string M_datafileSection;

    boost::shared_ptr<Epetra_FECrsMatrix> M_X;
    boost::shared_ptr<std::vector<long long int>> M_indexes_map_matrix_vectorized;
    boost::shared_ptr<Epetra_CrsMatrix> M_transposed_matrix;
    int M_offset;
};

} // namespace LifeV

#endif //  IMPORTMATRIXSNAPSHOTS_H_
