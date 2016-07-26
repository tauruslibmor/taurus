//@HEADER
/*!
    @file
    @author D. Forti
    @date 10-11-2015
*/

#ifndef IMPORTRHSSNAPSHOTS_H_
#define IMPORTRHSSNAPSHOTS_H_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

#include "EpetraExt_Utils.h"

namespace LifeV
{

class ImportRhsSnapshots
{
public:

	typedef Epetra_Comm comm_Type;
	typedef MapEpetra map_Type;
    typedef boost::shared_ptr<Epetra_Map> mapEpetraPtr_Type;
    typedef boost::shared_ptr<MapEpetra> mapPtr_Type;
	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
	typedef boost::shared_ptr<Epetra_FECrsMatrix> FECrsPtr_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;


	ImportRhsSnapshots ( const mapPtr_Type& map_rows, const commPtr_Type& communicator, const int numParameters, const GetPot& dataFile, const std::string& dataFile_section );

	ImportRhsSnapshots ( const int numParameters, const GetPot& dataFile);

    ~ImportRhsSnapshots() {};

    void readSnapshots();

    void buildColumnMap();

    void getMatrix ( FECrsPtr_Type& matrix );

    void buildRowMap ( );

private:

    int M_numParameters;
    GetPot M_datafile;
    mapEpetraPtr_Type M_map_column;
    mapPtr_Type M_map_rows;
    FECrsPtr_Type M_X;
    commPtr_Type M_comm;
    std::string M_datafileSection;

};

} // namespace LifeV

#endif //  IMPORTRHSSNAPSHOTS_H_
