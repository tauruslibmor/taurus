#ifndef DENSEHDF5_H
#define DENSEHDF5_H

#include "hdf5.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include <string>


class DenseHDF5
{
public:

	DenseHDF5();

	~DenseHDF5();

	bool Open(const std::string& filename);

	bool OpenParallel(const Epetra_Comm& Comm, const std::string& filename, int AccessType=H5F_ACC_RDWR);

	bool Create(const std::string& filename);

	void Close();

	bool IsOpen() const;

	bool WriteMatrixInt(const std::string& group, const std::string& dataset, int nx, int ny, const void * data);

	bool WriteMatrixDouble(const std::string& group, const std::string& dataset, int nx, int ny, const void * data);
    
    bool WriteEpetraSerialDenseMatrix(const std::string& group, const std::string& dataset, const Epetra_SerialDenseMatrix matrix);

	bool WriteVectorInt(const std::string& group, const std::string& dataset, int nx, const void * data);
    
    bool WriteVectorInt(const std::string& group, const std::string& dataset, const std::vector<int>& vec);
    
    bool WriteVectorLongInt(const std::string& group, const std::string& dataset, int nx, const void * data);
    
    bool WriteVectorLongInt(const std::string& group, const std::string& dataset, const std::vector<long long int>& vec);

	bool WriteVectorDouble(const std::string& group, const std::string& dataset, int nx, const void * data);
    
    bool WriteVectorDouble(const std::string& group, const std::string& dataset, const std::vector<double>& vec);

    bool WriteIntValue(const std::string& group, const std::string& dataset, const int v);

    bool OverWriteIntValue(const std::string& group, const std::string& dataset, const int v);

    bool WriteDoubleValue(const std::string& group, const std::string& dataset, const double v);

    bool ReadIntValue(const std::string& group, const std::string& dataset, int& v);
    
    bool ReadVectorProperties(const std::string& group, const std::string& dataset, void* nx );
    
    bool ReadMatrixProperties(const std::string& group, const std::string& dataset, int& nx, int& ny );

    bool ReadVectorLongInt(const std::string& group, const std::string& dataset, std::vector<long long int>& v);

    bool ReadVectorInt(const std::string& group, const std::string& dataset, std::vector<int>& v);
    
    bool ReadEpetraSerialDenseMatrix(const std::string& group, const std::string& dataset, Epetra_SerialDenseMatrix& matrix);


private:

    // methods
	hid_t CreateGroup(const std::string& group);

	bool WriteMatrixProperties(const std::string& group, const void * nx, const void * ny);

	bool WriteVectorProperties(const std::string& group, const void * nx);
    
    bool ReadInt(const std::string& group, const std::string& dataset, void* data);
    
    bool ReadLongInt(const std::string& group, const std::string& dataset, void* data);
    
    bool ReadDouble(const std::string& group, const std::string& dataset, void* data);

    // members
	bool M_isOpen;
	std::string M_filename;
	hid_t M_file_id;
	//const Epetra_Comm& M_comm;
};

#endif // DENSEHDF5_H
