//@HEADER
/*!
    @file
    @author D. Forti
    @date 10-11-2015
*/

#ifndef VECTORIZEDTOMATRICES_H_
#define VECTORIZEDTOMATRICES_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_Utils.h"
#include <EpetraExt_Transpose_RowMatrix.h>

namespace LifeV
{

class VectorizedToMatrices
{
public:

	typedef RegionMesh< LinearTetra > mesh_Type;
	typedef boost::shared_ptr<Epetra_Map> mapEpetraPtr_Type;

	VectorizedToMatrices ();

    ~VectorizedToMatrices() {};

    void setVectorizedMatrix ( const boost::shared_ptr<Epetra_FECrsMatrix>& vectorized_matrix) { M_vectorizedMatrix = vectorized_matrix; } ;

    void setMatrixRowMap ( const mapEpetraPtr_Type& map ) { M_map_matrix_rows = map; };

    void getFullMatrix ( const int& column_index, boost::shared_ptr<Epetra_FECrsMatrix>& matrix );

    void setOffset( const int& offset ) { M_offset = offset; };

private:

    void transpose();

    void reshapeIndex(const long long int& index_in, int& index_row, int& index_col );

    boost::shared_ptr<Epetra_FECrsMatrix> M_vectorizedMatrix;
    boost::shared_ptr<Epetra_CrsMatrix> M_vectorizedMatrix_transposed;
    mapEpetraPtr_Type M_map_matrix_rows;
    int M_offset;
    bool hasTranspose;

};

} // namespace LifeV

#endif //  VECTORIZEDTOMATRICES_H_
