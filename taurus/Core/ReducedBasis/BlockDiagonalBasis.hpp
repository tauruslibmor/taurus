//@HEADER
/*!
    @file
    @author F. Negri
    @date 27-01-2016
*/

#ifndef BDB_H_
#define BDB_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include <taurus/Core/Utilities/DenseHDF5.hpp>


namespace LifeV
{

class BlockDiagonalBasis
{
public:

	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;

	// constructor
	BlockDiagonalBasis ( const commPtr_Type& communicator,
						 const boost::shared_ptr<Epetra_FECrsMatrix>& V1,
						 const int num_rows_V1,
						 const boost::shared_ptr<Epetra_FECrsMatrix>& V2,
						 const boost::shared_ptr<Epetra_Map>& target_RowMap );

	// destructor
    ~BlockDiagonalBasis() {};

    void perform();

    void blockRow();

    void write( EpetraExt::HDF5& HDF5_exporter, const std::string& name);

    void getBasis ( boost::shared_ptr<Epetra_FECrsMatrix>& matrix );

    int getNumBases() const {return M_num_cols_V;};


private:

    commPtr_Type M_comm;
    boost::shared_ptr<Epetra_Map> M_map_rows;
    boost::shared_ptr<Epetra_FECrsMatrix> M_V1;
    boost::shared_ptr<Epetra_FECrsMatrix> M_V2;
    boost::shared_ptr<Epetra_FECrsMatrix> M_V;
    int M_num_cols_V;
    int M_num_rows_V1;
};

} // namespace Taurus

#endif //  REDUCTIONVECTOR_H_
