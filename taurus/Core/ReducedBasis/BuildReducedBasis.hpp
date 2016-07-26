//@HEADER
/*!
    @file
    @author F. Negri
    @date 27-01-2016
*/

#ifndef BUILDREDUCEDBASIS_H_
#define BUILDREDUCEDBASIS_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include <taurus/Core/Utilities/DenseHDF5.hpp>
#include <taurus/Core/ReducedBasis/V_POD.hpp>
#include <taurus/Core/ReducedBasis/ImportSnapshots.hpp>


namespace LifeV
{

class BuildReducedBasis
{
public:

	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<MatrixEpetra<Real>> matrixPtr_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;
	typedef RegionMesh< LinearTetra > mesh_Type;
	typedef FESpace<mesh_Type, MapEpetra>   fespace_Type;
	typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;
	typedef boost::shared_ptr<BCHandler> bcPtr_Type;

	// constructor
	BuildReducedBasis (const fespacePtr_Type& fe_space, bcPtr_Type bc, const GetPot& dataFile, const std::string& dataFile_section);

	// destructor
    ~BuildReducedBasis() {};

    // read snapshots, perform POD and DEIM
    void perform(const int NumSnapshots, const std::string solution_type="scalar");

    // export to HDF5 file
    void write(EpetraExt::HDF5& HDF5_exporter);
    void write(EpetraExt::HDF5& HDF5_exporter, DenseHDF5& HDF5dense_exporter);

    // getters
    int getNumBases() const {return M_NumPODbases;};
    int getNumSnapshots() const {return M_NumSnapshots;};
    std::vector<double> getSingularValues() const {return M_SingularValues;};
    void getPODbasis ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix );



private:

    commPtr_Type M_comm;
    GetPot M_datafile;
    fespacePtr_Type M_FESpace;
    boost::shared_ptr<BCHandler> M_bc_markingDOFs;

    std::string M_datafileSection;
    matrixPtr_Type M_Snapshots;
    boost::shared_ptr<Epetra_FECrsMatrix> M_PODbasis;
    int M_NumPODbases;
    int M_NumSnapshots;
    std::vector<double> M_SingularValues;
    std::string M_prefix;

};

} // namespace Taurus

#endif //  BUILDREDUCEDBASIS_H_
