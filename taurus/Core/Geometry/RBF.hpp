#ifndef RBF_H
#define RBF_H

#include "Epetra_MpiComm.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include <string>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

namespace LifeV
{

class RBF
{
public:

	typedef RegionMesh< LinearTetra > mesh_Type;
    typedef boost::shared_ptr<mesh_Type>  meshPtr_Type;
    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr< comm_Type > commPtr_Type;

	RBF();

	~RBF();

	void Setup(const std::string FileName);

	void Build(const std::string RBF_type, const std::vector<double>& ParamValue, const double Radius[]);

	void Build(const std::string RBF_type, const std::vector<double>& ParamValue, const double Radius[], const commPtr_Type& communicator);

	void MoveMesh(meshPtr_Type& mesh);

	void MoveMesh(meshPtr_Type& mesh, vectorPtr_Type& displacement );

	void ExportReferencePoints(std::string FilePrefix);

	void ExportDeformedPoints(std::string FilePrefix);

	void ResetMesh(meshPtr_Type& mesh, const meshPtr_Type& mesh_ref );

private:

    // methods
    double RBF_function(const double& distance, const double& radius);
    double distance(const int& i, const int& j);
    double distance(std::vector<double>& x, std::vector<double>& y);

    Epetra_SerialDenseMatrix Build_Matrix(const double &radius);

    // members
	std::string M_file;
	std::string M_RBFtype;
	double* M_ConstRadius;
	double M_NumPoints;
	double M_NumParam;
	std::vector<double> M_InterpValues;
	std::vector<std::vector<double> > M_InterpPoints;
	std::vector<std::vector<int> > M_ParamPoints;
	std::vector<std::vector<double> > M_Displacement_InterpPoints;
	Epetra_SerialDenseMatrix M_Matrix_x;
	Epetra_SerialDenseMatrix M_Matrix_y;
	Epetra_SerialDenseMatrix M_Matrix_z;
	Epetra_SerialDenseVector M_rhs_x;
	Epetra_SerialDenseVector M_rhs_y;
	Epetra_SerialDenseVector M_rhs_z;
	Epetra_SerialDenseVector M_weights_x;
	Epetra_SerialDenseVector M_weights_y;
	Epetra_SerialDenseVector M_weights_z;

};

} // namespace LifeV

#endif // RBF_H
