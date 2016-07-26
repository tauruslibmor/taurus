#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <lifev/core/LifeV.hpp>
#include "Epetra_Map.h"

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <taurus/Core/Utilities/DOF_Extractor.hpp>

namespace LifeV
{
    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    
    typedef MapEpetra map_Type;
    typedef boost::shared_ptr<map_Type> mapPtr_Type;
    
    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
    
class RTRed_Utils
{
public:
    
    RTRed_Utils(){}
    ~RTRed_Utils(){}
    
    void projectMatrix( const boost::shared_ptr<Epetra_FECrsMatrix>& A, const boost::shared_ptr<Epetra_FECrsMatrix>& V, const int num_bases,
                        const mapPtr_Type& map_rows, const mapPtr_Type& map_cols,  boost::shared_ptr<Epetra_FECrsMatrix>& out );
    
    void projectMatrix( const boost::shared_ptr<Epetra_FECrsMatrix>& A, const boost::shared_ptr<Epetra_FECrsMatrix>& V, const int num_bases,
                        boost::shared_ptr<Epetra_FECrsMatrix>& out );

    void transpose( boost::shared_ptr<Epetra_FECrsMatrix>& input_matrix, boost::shared_ptr<Epetra_CrsMatrix>& output_matrix);

    void extractRow(const boost::shared_ptr<Epetra_CrsMatrix>& matrix, const int row_index, double* srcValues, int* colIndices);

    void FECrsMatrixToMultiVector(const boost::shared_ptr<Epetra_FECrsMatrix>& input_matrix, boost::shared_ptr<Epetra_MultiVector>& vector, int numcols);

private:
    // matrixPtr_Type M_V;
    
}; // namespace RTRed_Utils
    
} // na

#endif //  UTILITIES_H_
