//@HEADER
/*!
    @file
    @author F. Negri
    @date 27-01-2016
*/

#ifndef HYPERREDUCTIONVECTOR_H_
#define HYPERREDUCTIONVECTOR_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include <taurus/Core/Utilities/DenseHDF5.hpp>


namespace LifeV
{

class HyperReductionVector
{
public:

	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<MatrixEpetra<Real>> matrixPtr_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;

	// constructor
	HyperReductionVector (const GetPot& dataFile, const std::string& dataFile_section, const commPtr_Type& communicator, const mapPtr_Type& map_rows);


	// destructor
    ~HyperReductionVector() {};

    // read snapshots, perform POD and DEIM
    void perform();

    // Galerkin projection onto a basis V
    void project(const boost::shared_ptr<Epetra_FECrsMatrix>& V, const int numColsV, EpetraExt::HDF5& HDF5_exporter);


    // export to HDF5 file
    void write(DenseHDF5& HDF5_exporter);

    // save Snapshots matrix to .m file
    void spySnapshots();

    // save matrix containing POD basis to .m file
    void spyPOD();

    // getters
    int getNumBases() const {return M_NumPODbases;};
    int getNumSnapshots() const {return M_NumSnapshots;};
    std::vector<double> getSingularValues() const {return M_SingularValues;};
    std::vector<int> getDeimGID() const {return M_deim_GID;};

    void buildMaskMatrix(boost::shared_ptr<Epetra_FECrsMatrix>& P);
    void DEIMjacobianLeftProjection(const Epetra_CrsMatrix* V, boost::shared_ptr<Epetra_FECrsMatrix>& matrix);
    void DEIMjacobianLeftProjection(const  boost::shared_ptr<Epetra_FECrsMatrix>& V, EpetraExt::HDF5& HDF5_exporter, std::string dataset="DEIMjacobianLeftProjection");


private:
    
    //void FECrs_to_SerialDense(const boost::shared_ptr<Epetra_FECrsMatrix>& A, Epetra_SerialDenseMatrix& A_dense );
    void SerialDense_to_FECrs(const Epetra_SerialDenseMatrix* A_dense, boost::shared_ptr<Epetra_FECrsMatrix>& A );

    commPtr_Type M_comm;
    GetPot M_datafile;
    std::string M_datafileSection;
    mapPtr_Type M_map_rows;
    boost::shared_ptr<Epetra_FECrsMatrix> M_Snapshots;
    boost::shared_ptr<Epetra_FECrsMatrix> M_PODbasis;
    int M_NumPODbases;
    int M_NumSnapshots;
    std::vector<double> M_SingularValues;
    std::vector<int> M_deim_GID;
    std::string M_HDF5outputFile;
    Epetra_SerialDenseMatrix M_DEIM_Phi;
    std::string M_prefix;
    
    boost::shared_ptr<Epetra_FECrsMatrix> M_maskMatrix;

};

} // namespace Taurus

#endif //  HYPERREDUCTIONVECTOR_H_
