//@HEADER
/*!
    @file
    @author F. Negri
    @date 27-01-2016
*/

#ifndef REDUCTIONVECTOR_H_
#define REDUCTIONVECTOR_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include <taurus/Core/Utilities/DenseHDF5.hpp>


namespace LifeV
{

class ReductionVector
{
public:

	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<MatrixEpetra<Real>> matrixPtr_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;

	// constructor
	ReductionVector (const GetPot& dataFile, const std::string& dataFile_section, const commPtr_Type& communicator, const mapPtr_Type& map_rows);


	// destructor
    ~ReductionVector() {};

    void setWeightMatrix ( boost::shared_ptr<Epetra_FECrsMatrix> matrix );

    // read snapshots, perform POD and DEIM
    void perform();

    // export to HDF5 file
    void write(EpetraExt::HDF5& HDF5_exporter, bool usePrefix=false);
    void write(EpetraExt::HDF5& HDF5_exporter, DenseHDF5& HDF5dense_exporter, bool usePrefix=false);

    // getters
    int getNumBases() const {return M_NumPODbases;};
    int getNumSnapshots() const {return M_NumSnapshots;};
    std::vector<double> getSingularValues() const {return M_SingularValues;};
    void getPODbasis ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix );


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
    std::string M_HDF5outputFile;
    std::string M_prefix;

    boost::shared_ptr<Epetra_FECrsMatrix> M_weightMatrix;
    bool M_useWeightMatrix;

};

} // namespace Taurus

#endif //  REDUCTIONVECTOR_H_
