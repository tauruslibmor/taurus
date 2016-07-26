//@HEADER
/*!
    @file
    @author F. Negri
    @date 27-01-2016
*/

#ifndef HYPERREDUCTIONMATRIX_H_
#define HYPERREDUCTIONMATRIX_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include <taurus/Core/Utilities/DenseHDF5.hpp>
#include "EpetraExt_HDF5.h"



namespace LifeV
{

class HyperReductionMatrix
{
public:

	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<MatrixEpetra<Real>> matrixPtr_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;

	// constructor
	HyperReductionMatrix (const GetPot& dataFile, const std::string& dataFile_section, const commPtr_Type& communicator, const mapPtr_Type& map_rows);


	// destructor
    ~HyperReductionMatrix() {};

    // read snapshots, perform POD and DEIM
    void perform();

    // export to HDF5 file
    void write(DenseHDF5& HDF5_exporter);

    // Galerkin projection onto a basis V
    void project(const boost::shared_ptr<Epetra_FECrsMatrix>& V, const int numColsV, EpetraExt::HDF5& HDF5_exporter);

    // save Snapshots matrix to .m file
    void spySnapshots();

    // save matrix containing POD basis to .m file
    void spyPOD();

    // getters
    void getPODbasis ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix );
    int getNumBases() const {return M_NumPODbases;};
    int getNumSnapshots() const {return M_NumSnapshots;};
    std::vector<double> getSingularValues() const {return M_SingularValues;};
    std::vector<long long int> getDeimGID() const {return M_deim_GID;};
    std::vector<int> getDeimDofs() const {return M_deim_dofs_GID;};
    std::vector<int> getDeimRows() const {return M_deim_row_indices;};
    std::vector<int> getDeimCols() const {return M_deim_col_indices;};
    int getOffset() const {return M_offset;};



private:

    commPtr_Type M_comm;
    GetPot M_datafile;
    std::string M_datafileSection;
    mapPtr_Type M_map_rows;
    boost::shared_ptr<Epetra_FECrsMatrix> M_Snapshots;
    boost::shared_ptr<Epetra_FECrsMatrix> M_PODbasis;
    int M_NumPODbases;
    int M_NumSnapshots;
    std::vector<double> M_SingularValues;
    std::vector<long long int> M_deim_GID;
    std::vector<int> M_deim_dofs_GID;
    std::vector<int> M_deim_row_indices;
    std::vector<int> M_deim_col_indices;
    std::string M_HDF5outputFile;
    Epetra_SerialDenseMatrix M_DEIM_Phi;
    std::string M_prefix;
    int M_offset;

};

} // namespace Taurus

#endif //  HYPERREDUCTIONMATRIX_H_
