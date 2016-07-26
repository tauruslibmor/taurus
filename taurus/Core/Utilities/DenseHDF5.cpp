#include <taurus/Core/Utilities/DenseHDF5.hpp>
#include "EpetraExt_Exception.h"
#include "Epetra_MpiComm.h"

#define CHECK_STATUS(status) \
      { if (status < 0) \
        throw(EpetraExt::Exception(__FILE__, __LINE__, \
                        "DenseHDF5: function H5Giterater returned a negative value")); }

using namespace std;

//=========================================================================
// Constructor //
DenseHDF5::DenseHDF5() :
				M_file_id  ( -1 ),
				M_isOpen   ( false )
{}
//=========================================================================

// Destructor //
DenseHDF5::~DenseHDF5()
{
	Close();
}
//=========================================================================
bool DenseHDF5::Create(const string& filename)
{
	if (M_file_id >= 0)
		Close();

	M_filename = filename;

	M_file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	return M_file_id >= 0;
}
//=========================================================================
bool DenseHDF5::Open(const string& filename)
{
	if (M_file_id >= 0)
		Close();

	M_filename = filename;

	//plist_id_ = H5Pcreate(H5P_FILE_ACCESS);
	//M_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id_);

	M_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

	return M_file_id >= 0;
}
//=========================================================================
bool DenseHDF5::OpenParallel(const Epetra_Comm& Comm, const string& filename, int AccessType)
{

	//M_comm = Comm;

	if (IsOpen())
		throw(EpetraExt::Exception(__FILE__, __LINE__,
				"an HDF5 is already open, first close the current one",
				"using method Close(), then open/create a new one"));

	M_filename = filename;

	// Set up file access property list with parallel I/O access
	hid_t plist_id_ = H5Pcreate(H5P_FILE_ACCESS);

#ifdef HAVE_MPI
	// Create property list for collective dataset write.
	const Epetra_MpiComm* MpiComm ( dynamic_cast<const Epetra_MpiComm*> (&Comm) );

	if (MpiComm == 0)
		H5Pset_fapl_mpio(plist_id_, MPI_COMM_WORLD, MPI_INFO_NULL);
	else
		H5Pset_fapl_mpio(plist_id_, MpiComm->Comm(), MPI_INFO_NULL);
#endif

// create the file collectively and release property list identifier.
	M_file_id = H5Fopen(filename.c_str(), AccessType, plist_id_);
	H5Pclose(plist_id_);

	return M_file_id >= 0;
}
//=========================================================================
void DenseHDF5::Close()
{
	if (M_file_id >= 0)
	{
		H5Fclose(M_file_id);
		M_file_id = -1;
	}
}
//=========================================================================
bool DenseHDF5::IsOpen() const
{
	return M_file_id >= 0;
}

hid_t DenseHDF5::CreateGroup(const string& group)
{
	hid_t group_id = -1;

	if ( group != "/")
	{
		if ( H5Lexists(M_file_id, group.c_str(), H5P_DEFAULT)==0 )
		{
			group_id = H5Gcreate(M_file_id, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		else
		{
			group_id = H5Gopen(M_file_id, group.c_str(), H5P_DEFAULT);
		}
	}

return group_id;
}
//=========================================================================
bool DenseHDF5::WriteVectorProperties(const std::string& group, const void * nx)
{

	if (M_file_id < 0)
		return false;
	// open group
	hid_t group_id = CreateGroup(group);

	// Create the data space for the dataset
	hsize_t     dims[1];
	dims[0]     = 1;

	hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

	// Create the dataset
	string full_dset_name = group+"/VecLength";
	hid_t dataset_id = H5Dcreate(M_file_id, full_dset_name.c_str(), H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write the dataset
	herr_t      status;

	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nx);
	CHECK_STATUS(status);

	// End access to the dataset and release resources used by it
	status = H5Dclose(dataset_id);
	CHECK_STATUS(status);

	// Terminate access to the data space
	status = H5Sclose(dataspace_id);
	CHECK_STATUS(status);

	return status >= 0;
}
//=========================================================================
bool DenseHDF5::WriteMatrixProperties(const std::string& group, const void * nx, const void * ny)
{
    
    if (M_file_id < 0)
    return false;
    // open group
    hid_t group_id = CreateGroup(group);
    
    // Create the data space for the dataset
    hsize_t     dims[1];
    dims[0]     = 1;
    
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    
    // Create the dataset
    string full_dset_name = group+"/RowNum";
    hid_t dataset_id = H5Dcreate(M_file_id, full_dset_name.c_str(), H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write the dataset
    herr_t      status;
    
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nx);
    CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    // Terminate access to the data space
    status = H5Sclose(dataspace_id);
    CHECK_STATUS(status);
    
    // Create the data space for the dataset
    dataspace_id = H5Screate_simple(1, dims, NULL);
    
    // Create the dataset
    full_dset_name = group+"/ColNum";
    dataset_id = H5Dcreate(M_file_id, full_dset_name.c_str(), H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ny);
    CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    // Terminate access to the data space
    status = H5Sclose(dataspace_id);
    CHECK_STATUS(status);
    
    return status >= 0;
}
//=========================================================================
bool DenseHDF5::WriteMatrixInt(const string& group, const string& dataset, int nx, int ny, const void * data)
{
	if (M_file_id < 0)
		return false;

    // open group
	hid_t group_id = CreateGroup(group);
    string full_dset_name = group + dataset;
    group_id = CreateGroup(full_dset_name);

    // WriteMatrixProperties
    WriteMatrixProperties(full_dset_name, &nx, &ny);
    
    if (nx == 0 && ny == 0)
    	return true;// ragionarla sulla scrittura di strutture vuote, puo servire?
    if (nx < 0 || ny < 0)
        return false;// ragionarla sulla scrittura di strutture vuote, puo servire?
    
	// Create the data space for the dataset
	hsize_t     dims[2];
	dims[0] = nx;
	dims[1] = ny;

	hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

	// Create the dataset
	string full_dset_name2 = full_dset_name+"/Values";
	hid_t dataset_id = H5Dcreate(M_file_id, full_dset_name2.c_str(), H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write the dataset
	herr_t      status;

	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	CHECK_STATUS(status);

	// End access to the dataset and release resources used by it
	status = H5Dclose(dataset_id);
	CHECK_STATUS(status);

	// Terminate access to the data space
	status = H5Sclose(dataspace_id);
	CHECK_STATUS(status);

	return status >= 0;
}
//=========================================================================
bool DenseHDF5::WriteMatrixDouble(const string& group, const string& dataset, int nx, int ny, const void * data)
{
    if (M_file_id < 0)
    return false;
    
    // open group
    hid_t group_id = CreateGroup(group);
    string full_dset_name = group + dataset;
    group_id = CreateGroup(full_dset_name);
    
    // WriteMatrixProperties
    WriteMatrixProperties(full_dset_name, &nx, &ny);
    
    if (nx == 0 && ny == 0)
    return true;// ragionarla sulla scrittura di strutture vuote, puo servire?
    if (nx < 0 || ny < 0)
    return false;// ragionarla sulla scrittura di strutture vuote, puo servire?
    
    // Create the data space for the dataset
    hsize_t     dims[2];
    dims[0] = nx;
    dims[1] = ny;
    
    hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
    
    // Create the dataset
    string full_dset_name2 = full_dset_name+"/Values";
    hid_t dataset_id = H5Dcreate(M_file_id, full_dset_name2.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write the dataset
    herr_t      status;
    
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    // Terminate access to the data space
    status = H5Sclose(dataspace_id);
    CHECK_STATUS(status);
    
    return status >= 0;
}
//=========================================================================
bool DenseHDF5::WriteEpetraSerialDenseMatrix(const std::string& group, const std::string& dataset, const Epetra_SerialDenseMatrix matrix)
{
 
    double data[matrix.M()][matrix.N()];
    for (int i = 0; i < matrix.M(); i++)
       for (int j = 0; j < matrix.N(); j++)
           data[i][j] = matrix(i,j);
    
    return WriteMatrixDouble(group, dataset, matrix.M(), matrix.N(), data);
    
}
//=========================================================================
bool DenseHDF5::WriteVectorInt(const string& group, const string& dataset, int nx, const void * data)
{
	if (M_file_id < 0)
		return false;

	// open group
	hid_t group_id = CreateGroup(group);

	// open group
	string full_dset_name = group + dataset;
	group_id = CreateGroup(full_dset_name);

	// WriteVectorProperties
	WriteVectorProperties(full_dset_name, &nx);

	// Create the data space for the dataset
	hsize_t     dims[1];
	dims[0] = nx;

	hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

	// Create the dataset
    string full_dset_name2 = full_dset_name+"/Values";
	hid_t dataset_id = H5Dcreate(M_file_id, full_dset_name2.c_str(), H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if (nx == 0)
		return true;// ragionarla sulla scrittura di strutture vuote, puo servire?

	if (nx < 0)
		return false;

	// Write the dataset
	herr_t   status;

	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	CHECK_STATUS(status);

	// End access to the dataset and release resources used by it
	status = H5Dclose(dataset_id);
	CHECK_STATUS(status);

	// Terminate access to the data space
	status = H5Sclose(dataspace_id);
	CHECK_STATUS(status);

	return status >= 0;
}
//=========================================================================
bool DenseHDF5::WriteVectorInt(const std::string& group, const std::string& dataset, const std::vector<int>& vec)
{
    int data[vec.size()];
    for (int i = 0; i < vec.size(); i++)
        data[i] = vec[i];
    
    return WriteVectorInt(group, dataset, (int)vec.size(), data);
}
//=========================================================================
bool DenseHDF5::WriteVectorDouble(const string& group, const string& dataset, int nx, const void * data)
{
    if (M_file_id < 0)
    return false;
    
    // open group
    hid_t group_id = CreateGroup(group);
    
    // open group
    string full_dset_name = group + dataset;
    group_id = CreateGroup(full_dset_name);
    
    // WriteVectorProperties
    WriteVectorProperties(full_dset_name, &nx);
    
    // Create the data space for the dataset
    hsize_t     dims[1];
    dims[0] = nx;
    
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    
    // Create the dataset
    string full_dset_name2 = full_dset_name+"/Values";
    hid_t dataset_id = H5Dcreate(M_file_id, full_dset_name2.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    if (nx == 0)
    return true;// ragionarla sulla scrittura di strutture vuote, puo servire?
    
    if (nx < 0)
    return false;
    
    // Write the dataset
    herr_t   status;
    
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    // Terminate access to the data space
    status = H5Sclose(dataspace_id);
    CHECK_STATUS(status);
    
    return status >= 0;
}
//=========================================================================
bool DenseHDF5::WriteVectorDouble(const std::string& group, const std::string& dataset, const std::vector<double>& vec)
{
    double data[vec.size()];
    for (int i = 0; i < vec.size(); i++)
    data[i] = vec[i];
    
    return WriteVectorDouble(group, dataset, (int)vec.size(), data);
}
//=========================================================================
bool DenseHDF5::WriteVectorLongInt(const string& group, const string& dataset, int nx, const void * data)
{
    if (M_file_id < 0)
    return false;
    
    // open group
    hid_t group_id = CreateGroup(group);
    
    // open group
    string full_dset_name = group + dataset;
    group_id = CreateGroup(full_dset_name);
    
    // WriteVectorProperties
    WriteVectorProperties(full_dset_name, &nx);
    
    // Create the data space for the dataset
    hsize_t     dims[1];
    dims[0] = nx;
    
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    
    // Create the dataset
    string full_dset_name2 = full_dset_name+"/Values";
    hid_t dataset_id = H5Dcreate(M_file_id, full_dset_name2.c_str(), H5T_NATIVE_LLONG, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    if (nx == 0)
    return true;// ragionarla sulla scrittura di strutture vuote, puo servire?
    
    if (nx < 0)
    return false;
    
    // Write the dataset
    herr_t   status;
    
    status = H5Dwrite(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    // Terminate access to the data space
    status = H5Sclose(dataspace_id);
    CHECK_STATUS(status);
    
    return status >= 0;
}
//=========================================================================
bool DenseHDF5::WriteVectorLongInt(const std::string& group, const std::string& dataset, const std::vector<long long int>& vec)
{
    long long int data[vec.size()];
    for (int i = 0; i < vec.size(); i++)
         data[i] = vec[i];
    
    return WriteVectorLongInt(group, dataset, (int)vec.size(), data);
}
//=========================================================================
bool DenseHDF5::WriteIntValue(const string& group, const string& dataset, const int v)
{
    if (M_file_id < 0)
    return false;
    
    // open group
    hid_t group_id = CreateGroup(group);
    
    // Create the data space for the dataset
    hsize_t     dims[1];
    dims[0]     =  1;
    
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    
    // Create the dataset
    string full_dset_name = group + dataset;
    hid_t dataset_id = H5Dcreate(M_file_id, full_dset_name.c_str(), H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write the dataset
    herr_t   status;
    
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
    CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    // Terminate access to the data space
    status = H5Sclose(dataspace_id);
    CHECK_STATUS(status);
    
    return status >= 0;
}
//=========================================================================
bool DenseHDF5::OverWriteIntValue(const string& group, const string& dataset, const int v)
{
    if (M_file_id < 0)
    return false;

    // open group
    hid_t group_id = CreateGroup(group);

    // Create the data space for the dataset
    hsize_t     dims[1];
    dims[0]     =  1;

    // Create the dataset
    string full_dset_name = group + dataset;
    //hid_t plist = H5Pcreate (H5P_DATATYPE_ACCESS);
    hid_t dataset_id = H5Dopen(M_file_id, full_dset_name.c_str(), H5P_DEFAULT);

    // Write the dataset
    herr_t   status;

    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
    CHECK_STATUS(status);

    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);

    // Terminate access to the data space
    CHECK_STATUS(status);

    return status >= 0;
}
//=========================================================================
bool DenseHDF5::WriteDoubleValue(const string& group, const string& dataset, const double v)
{
    if (M_file_id < 0)
    return false;
    
    // open group
    hid_t group_id = CreateGroup(group);
    
    // Create the data space for the dataset
    hsize_t     dims[1];
    dims[0]     =  1;
    
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    
    // Create the dataset
    string full_dset_name = group + dataset;
    hid_t dataset_id = H5Dcreate(M_file_id, full_dset_name.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write the dataset
    herr_t   status;
    
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
    CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    // Terminate access to the data space
    status = H5Sclose(dataspace_id);
    CHECK_STATUS(status);
    
    return status >= 0;
}
//=========================================================================
bool DenseHDF5::ReadIntValue(const string& group, const string& dataset, int& v)
{
    if (M_file_id < 0)
    return false;

    herr_t   status;
    
    hid_t group_id   = H5Gopen(M_file_id, group.c_str(), H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(group_id, dataset.c_str(), H5P_DEFAULT);
    
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
    
    CHECK_STATUS(status);
    
    return status >= 0;
}
//=========================================================================
bool DenseHDF5::ReadInt(const string& group, const string& dataset, void* data)
{
    if (M_file_id < 0)
    return false;

	herr_t   status;

	hid_t group_id = H5Gopen(M_file_id, group.c_str(), H5P_DEFAULT);
	hid_t dataset_id = H5Dopen(group_id, dataset.c_str(), H5P_DEFAULT);

	status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);

	return status >= 0;
}
//=========================================================================
bool DenseHDF5::ReadLongInt(const string& group, const string& dataset, void* data)
{
    if (M_file_id < 0)
    return false;

    herr_t   status;
    
    hid_t group_id = H5Gopen(M_file_id, group.c_str(), H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(group_id, dataset.c_str(), H5P_DEFAULT);
    
    status = H5Dread(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    
    CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    return status >= 0;
}
//=========================================================================
bool DenseHDF5::ReadDouble(const string& group, const string& dataset, void* data)
{
    if (M_file_id < 0)
    return false;

	herr_t   status, status_n;

	hid_t group_id = H5Gopen(M_file_id, group.c_str(), H5P_DEFAULT);
	hid_t dataset_id = H5Dopen(group_id, dataset.c_str(), H5P_DEFAULT);
    
    /*
    ///////////////////////////////////////
    // Get dataset rank and dimension.
    hsize_t     dims[2];
    hsize_t     chunk_dims[2];

    hid_t filespace = H5Dget_space(dataset_id);
    int rank        = H5Sget_simple_extent_ndims(filespace);
    status_n        = H5Sget_simple_extent_dims(filespace, dims, NULL);
    printf("dataset rank %d, dimensions %lu x %lu\n",
           rank, (unsigned long)(dims[0]), (unsigned long)(dims[1]));
    
    // Get creation properties list.
    hid_t cparms = H5Dget_create_plist(dataset_id);
    
    // Check if dataset is chunked.
    if (H5D_CHUNKED == H5Pget_layout(cparms))  {
        
        // Get chunking information: rank and dimensions
        int rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
        printf("chunk rank %d, dimensions %lu x %lu\n", rank_chunk,
               (unsigned long)(chunk_dims[0]), (unsigned long)(chunk_dims[1]));
    }
    
    //Define the memory space to read dataset.
    hid_t memspace = H5Screate_simple(2,dims,NULL);
    
    // Read dataset back and display.
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                     H5P_DEFAULT, data);
    //printf("\n");
    //printf("Dataset: \n");
    //for (int j = 0; j < dims[0]; j++) {
    //    for (int i = 0; i < dims[1]; i++) printf("%d ", data[j][i]);
    //    printf("\n");
    //}
    //////////////////////////////////////
    */

	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);

	return status >= 0;
}
//=========================================================================
bool DenseHDF5::ReadVectorProperties(const string& group, const string& dataset, void* nx )
{
    if (M_file_id < 0)
    return false;

    herr_t status;
    
    string full_group_name = group + dataset;
    
    hid_t group_id = H5Gopen(M_file_id, full_group_name.c_str(), H5P_DEFAULT);
    
    string full_dset_name = full_group_name + "/VecLength";

    hid_t dataset_id = H5Dopen(group_id, full_dset_name.c_str(), H5P_DEFAULT);
    
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nx);
    
    CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    return status >= 0;
}
//=========================================================================
bool DenseHDF5::ReadVectorLongInt(const std::string& group, const std::string& dataset, std::vector<long long int>& v)
{
    if (M_file_id < 0)
    return false;
    
    int nx;
    ReadVectorProperties(group, dataset, &nx );
    
    long long int data[nx];
    
    ReadLongInt(group+dataset+"/", "Values", data);
    
    for (int i = 0; i < nx; i++)
        v.push_back( data[i] );
    return true;
    
}
//=========================================================================
bool DenseHDF5::ReadVectorInt(const std::string& group, const std::string& dataset, std::vector<int>& v)
{
    if (M_file_id < 0)
    return false;
    
    int nx;
    ReadVectorProperties(group, dataset, &nx );
    
    int data[nx];
    
    ReadInt(group+dataset+"/", "Values", data);
    
    for (int i = 0; i < nx; i++)
       v.push_back( data[i] );
    
    return true;
    
}
//=========================================================================
bool DenseHDF5::ReadMatrixProperties(const string& group, const string& dataset, int& nx, int& ny )
{
    if (M_file_id < 0)
    return false;
    
    herr_t status;
    
    string full_group_name = group + dataset;
    
    hid_t group_id = H5Gopen(M_file_id, full_group_name.c_str(), H5P_DEFAULT);
    
    string full_dset_name = full_group_name + "/RowNum";
    
    hid_t dataset_id = H5Dopen(group_id, full_dset_name.c_str(), H5P_DEFAULT);
    
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nx);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    full_dset_name = full_group_name + "/ColNum";
    
    dataset_id = H5Dopen(group_id, full_dset_name.c_str(), H5P_DEFAULT);
    
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ny);
    
    CHECK_STATUS(status);
    
    // End access to the dataset and release resources used by it
    status = H5Dclose(dataset_id);
    CHECK_STATUS(status);
    
    return status >= 0;
}
//=========================================================================
bool DenseHDF5::ReadEpetraSerialDenseMatrix(const string& group, const string& dataset, Epetra_SerialDenseMatrix& matrix)
{
    if (M_file_id < 0)
    return false;
    
    int nx, ny;
    ReadMatrixProperties(group, dataset, nx, ny );
    
    double data[nx][ny];
        
    matrix.Shape(nx, ny);
    
    ReadDouble(group+dataset+"/", "Values", data);
    
    for (int i = 0; i < nx; i++)
       for (int j = 0; j < ny; j++)
         matrix(i,j) = data[i][j];
 
    return true;
}
//=========================================================================

