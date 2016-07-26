#include <taurus/Core/Geometry/RBF.hpp>
#include <chrono>
using namespace std::chrono;

using namespace std;

using namespace LifeV;

//=========================================================================
// Constructor //
RBF::RBF() :
	M_NumPoints (0),
	M_NumParam  (0)
{
	M_ConstRadius = new double[3];
}
//=========================================================================
// Destructor //
RBF::~RBF()
{}
//=========================================================================
void RBF::Setup(const std::string FileName)
{
	M_file = FileName;

	std::fstream myfile(FileName, std::ios_base::in);

	myfile >> M_NumPoints;
	//cout << "M_NumPoints = " << M_NumPoints << endl;

	vector<vector<double> > InterpPoints(M_NumPoints, std::vector<double>(3));

	for (int i = 0; i < M_NumPoints; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			myfile >> InterpPoints[i][j];
			//printf("%f \n", InterpPoints[i][j]);
		}
	}

	M_InterpPoints = InterpPoints;

	myfile >> M_NumParam;
	//cout << "M_NumParam = " << M_NumParam << endl;

	vector<vector<int> > InterpParam(M_NumParam, std::vector<int>(2));

	for (int i = 0; i < M_NumParam; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			double tmp;
			myfile >> tmp;
			InterpParam[i][j] = (int)(tmp - 1.0);
			//printf("%d \n", InterpParam[i][j]);
		}
	}
	M_ParamPoints = InterpParam;

}
//=========================================================================
void RBF::Build(const std::string RBF_type, const std::vector<double>& ParamValue, const double Radius[])
{
	for (int i = 0; i < 3; i++)
		M_ConstRadius[i] = Radius[i];

	M_RBFtype = RBF_type;

	// create target displacements
	vector<vector<double> > Displacement_InterpPoints(M_NumPoints, std::vector<double>(3));

	for (int i = 0; i < M_NumPoints; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Displacement_InterpPoints[i][j] = 0;
		}
	}

	for (int i = 0; i < M_NumParam; i++)
	{
		Displacement_InterpPoints[ M_ParamPoints[i][0] ][ M_ParamPoints[i][1] ] = ParamValue[i];
	}

	M_Displacement_InterpPoints = Displacement_InterpPoints;

	// build interpolation matrix
	M_Matrix_x = Build_Matrix(Radius[0]);
	M_Matrix_y = Build_Matrix(Radius[1]);
	M_Matrix_z = Build_Matrix(Radius[2]);

	// build interpolation rhs
	M_rhs_x.Size(M_NumPoints+4);
	for (int i = 0; i < M_NumPoints; i++)
	{
		M_rhs_x(i) = Displacement_InterpPoints[i][0];
	}


	M_rhs_y.Size(M_NumPoints+4);
	for (int i = 0; i < M_NumPoints; i++)
	{
		M_rhs_y(i) = Displacement_InterpPoints[i][1];
	}

	M_rhs_z.Size(M_NumPoints+4);
	for (int i = 0; i < M_NumPoints; i++)
	{
		M_rhs_z(i) = Displacement_InterpPoints[i][2];
	}

	// compute coefficients
 	M_weights_x.Size( M_NumPoints+4 );
	Epetra_SerialDenseSolver my_solver;

	my_solver.SetMatrix( M_Matrix_x );
	my_solver.SetVectors( M_weights_x , M_rhs_x );
	my_solver.Solve();

	M_weights_y.Size( M_NumPoints+4 );
	my_solver.SetMatrix( M_Matrix_y );
	my_solver.SetVectors( M_weights_y , M_rhs_y );
	my_solver.Solve();

	M_weights_z.Size( M_NumPoints+4 );
	my_solver.SetMatrix( M_Matrix_z );
	my_solver.SetVectors( M_weights_z , M_rhs_z );
	my_solver.Solve();

}
//=========================================================================
void RBF::Build(const std::string RBF_type, const std::vector<double>& ParamValue, const double Radius[], const commPtr_Type& Comm)
{

	high_resolution_clock::time_point tStart;
	high_resolution_clock::time_point tEnd;


	for (int i = 0; i < 3; i++)
		M_ConstRadius[i] = Radius[i];

	M_RBFtype = RBF_type;

	// create target displacements
	vector<vector<double> > Displacement_InterpPoints(M_NumPoints, std::vector<double>(3));

	for (int i = 0; i < M_NumPoints; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Displacement_InterpPoints[i][j] = 0;
		}
	}

	for (int i = 0; i < M_NumParam; i++)
	{
		Displacement_InterpPoints[ M_ParamPoints[i][0] ][ M_ParamPoints[i][1] ] = ParamValue[i];
	}

	M_Displacement_InterpPoints = Displacement_InterpPoints;

	// build interpolation matrix
	M_Matrix_x = Build_Matrix(Radius[0]);
	M_Matrix_y = Build_Matrix(Radius[1]);
	M_Matrix_z = Build_Matrix(Radius[2]);

	// build interpolation rhs
	M_rhs_x.Size(M_NumPoints+4);
	for (int i = 0; i < M_NumPoints; i++)
	{
		M_rhs_x(i) = Displacement_InterpPoints[i][0];
	}


	M_rhs_y.Size(M_NumPoints+4);
	for (int i = 0; i < M_NumPoints; i++)
	{
		M_rhs_y(i) = Displacement_InterpPoints[i][1];
	}

	M_rhs_z.Size(M_NumPoints+4);
	for (int i = 0; i < M_NumPoints; i++)
	{
		M_rhs_z(i) = Displacement_InterpPoints[i][2];
	}

	// compute coefficients
	Epetra_SerialDenseVector root_weights_x;
	root_weights_x.Size( M_NumPoints+4 );

	Epetra_SerialDenseVector root_weights_y;
	root_weights_y.Size( M_NumPoints+4 );

	Epetra_SerialDenseVector root_weights_z;
	root_weights_z.Size( M_NumPoints+4 );

	int root = 0;

	if ( Comm->MyPID() == 0 )
	{
		Epetra_SerialDenseSolver my_solver;

		tStart = high_resolution_clock::now();

		my_solver.SetMatrix( M_Matrix_x );
		my_solver.SetVectors( root_weights_x , M_rhs_x );

		tEnd = high_resolution_clock::now();
					if ( Comm->MyPID() == 0 )
					std::cout << "\n\nElapsed time for compute coeffients X: " << duration_cast<milliseconds>( tEnd - tStart ).count() << " milliseconds\n\n";


		my_solver.Solve();



		my_solver.SetMatrix( M_Matrix_y );
		my_solver.SetVectors( root_weights_y , M_rhs_y );
		my_solver.Solve();

		my_solver.SetMatrix( M_Matrix_z );
		my_solver.SetVectors( root_weights_z , M_rhs_z );
		my_solver.Solve();
	}

	// broadcast

	Comm->Broadcast	(root_weights_x.Values(), M_NumPoints+4, root );
	Comm->Broadcast	(root_weights_y.Values(), M_NumPoints+4, root );
	Comm->Broadcast	(root_weights_z.Values(), M_NumPoints+4, root );

	M_weights_x.Size( M_NumPoints+4 );
	M_weights_y.Size( M_NumPoints+4 );
	M_weights_z.Size( M_NumPoints+4 );

	for (int i = 0; i < M_NumPoints; i++)
	{
		M_weights_x(i) = root_weights_x(i);
		M_weights_y(i) = root_weights_y(i);
		M_weights_z(i) = root_weights_z(i);
	}
}
//=========================================================================
double RBF::RBF_function(const double& d, const double& r)
{

	double val = -1;

	if ( M_RBFtype.compare("gaussian")==0 )
	{
		val = exp(-0.5*d*d/(r*r));
	}

	if ( M_RBFtype.compare("thinplate")==0 )
	{
		val = d*d*log(d+1);
	}

	if ( M_RBFtype.compare("cubic")==0 )
	{
		val = (d*d*d);
	}

	return val;

}
//=========================================================================
double RBF::distance(const int& i, const int& j)
{

	double tmp = 0;

	for (int k = 0; k < 3; k++)
	{
		double tmp2 = (M_InterpPoints[i][k] - M_InterpPoints[j][k]);
		tmp += (tmp2*tmp2);
	}

	double d = sqrt(tmp);

	return d;
}
//=========================================================================
double RBF::distance(vector<double>& x, vector<double>& y)
{

	double tmp = 0;

	for (int k = 0; k < 3; k++)
	{
		double tmp2 = (x[k] - y[k]);
		tmp += (tmp2*tmp2);
	}

	double d = sqrt(tmp);

	return d;
}
//=========================================================================
Epetra_SerialDenseMatrix RBF::Build_Matrix(const double &radius)
{
	Epetra_SerialDenseMatrix A;
	A.Shape(M_NumPoints+4, M_NumPoints+4);


	for (int i = 0; i < M_NumPoints; i++)
	{
		for (int j = 0; j < M_NumPoints; j++)
		{
			double d  = distance(i,j);
			A(i,j)    = RBF_function(d, radius);
		}
	}


	for (int i = M_NumPoints; i < M_NumPoints+1; i++)
	{
		for (int j = 0; j < M_NumPoints; j++)
		{
			A(i,j)    = 1;
			A(j,i)    = A(i,j);
		}
	}


	int dim = 0;
	for (int i = M_NumPoints+1; i < M_NumPoints+4; i++)
	{
		for (int j = 0; j < M_NumPoints; j++)
		{
			A(i,j)    = M_InterpPoints[j][dim];
			A(j,i)    = A(i,j);
		}
		dim += 1;
	}


	return A;
}
//=========================================================================
void RBF::MoveMesh(meshPtr_Type& mesh)
{

	for (int i = 0; i < mesh->numVertices(); i++)
	{
		vector<double> this_point = {mesh->point(i).x(), mesh->point(i).y(), mesh->point(i).z()};
		//cout << "this_point = " << mesh->point(i).x() << ", " << mesh->point(i).y() << ", " << mesh->point(i).z() << "\n";

		// x-component
		double disp = 0;

		for (int j = 0; j < M_NumPoints; j++)
		{
			double d  = distance(M_InterpPoints[j], this_point);
			double tmp  = RBF_function(d, M_ConstRadius[0]);
			disp += M_weights_x(j)  * tmp;
		}

		disp +=  M_weights_x(M_NumPoints) ;

		for (int j = M_NumPoints+1; j < M_NumPoints+4; j++)
		{
			disp += M_weights_x(j) * this_point[j-M_NumPoints-1];
		}

		//cout << "i = " << i << ", disp = " << disp << "\n";
		mesh->point(i).x() += disp;

		// y-component
		disp = 0;

		for (int j = 0; j < M_NumPoints; j++)
		{
			double d  = distance(M_InterpPoints[j], this_point);
			double tmp  = RBF_function(d, M_ConstRadius[1]);
			disp += M_weights_y(j)  * tmp;
		}

		disp +=  M_weights_y(M_NumPoints) ;

		for (int j = M_NumPoints+1; j < M_NumPoints+4; j++)
		{
			disp += M_weights_y(j) * this_point[j-M_NumPoints-1];
		}

		mesh->point(i).y() += disp;

		// z-component
		disp = 0;

		for (int j = 0; j < M_NumPoints; j++)
		{
			double d  = distance(M_InterpPoints[j], this_point);
			double tmp  = RBF_function(d, M_ConstRadius[2]);
			disp += M_weights_z(j)  * tmp;
		}

		disp +=  M_weights_z(M_NumPoints) ;

		for (int j = M_NumPoints+1; j < M_NumPoints+4; j++)
		{
			disp += M_weights_z(j) * this_point[j-M_NumPoints-1];
		}

		mesh->point(i).z() += disp;

	}

}
//=========================================================================

void RBF::MoveMesh(meshPtr_Type& mesh, vectorPtr_Type& displacement )
{
	for (int i = 0; i < mesh->numVertices(); i++)
	{
		int id = mesh->point(i).id();

		vector<double> this_point = {mesh->point(i).x(), mesh->point(i).y(), mesh->point(i).z()};
		//cout << "this_point = " << mesh->point(i).x() << ", " << mesh->point(i).y() << ", " << mesh->point(i).z() << "\n";

		// x-component
		double disp = 0;

		for (int j = 0; j < M_NumPoints; j++)
		{
			double d  = distance(M_InterpPoints[j], this_point);
			double tmp  = RBF_function(d, M_ConstRadius[0]);
			disp += M_weights_x(j)  * tmp;
		}

		disp +=  M_weights_x(M_NumPoints) ;

		for (int j = M_NumPoints+1; j < M_NumPoints+4; j++)
		{
			disp += M_weights_x(j) * this_point[j-M_NumPoints-1];
		}

		//cout << "i = " << i << ", disp = " << disp << "\n";

		if ( displacement->blockMap().LID ( id ) != -1 )
			(*displacement)[id] += disp;

		mesh->point(i).x() += disp;

		// y-component
		disp = 0;

		for (int j = 0; j < M_NumPoints; j++)
		{
			double d  = distance(M_InterpPoints[j], this_point);
			double tmp  = RBF_function(d, M_ConstRadius[1]);
			disp += M_weights_y(j)  * tmp;
		}

		disp +=  M_weights_y(M_NumPoints) ;

		for (int j = M_NumPoints+1; j < M_NumPoints+4; j++)
		{
			disp += M_weights_y(j) * this_point[j-M_NumPoints-1];
		}

		if ( displacement->blockMap().LID ( id + displacement->size()/3 ) != -1 )
			(*displacement)[ id +displacement->size()/3] += disp;

		mesh->point(i).y() += disp;

		// z-component
		disp = 0;

		for (int j = 0; j < M_NumPoints; j++)
		{
			double d  = distance(M_InterpPoints[j], this_point);
			double tmp  = RBF_function(d, M_ConstRadius[2]);
			disp += M_weights_z(j)  * tmp;
		}
		disp +=  M_weights_z(M_NumPoints) ;

		for (int j = M_NumPoints+1; j < M_NumPoints+4; j++)
		{
			disp += M_weights_z(j) * this_point[j-M_NumPoints-1];
		}

		if ( displacement->blockMap().LID ( id +2*displacement->size()/3 ) != -1 )
			(*displacement)[id+2*displacement->size()/3] += disp;

		mesh->point(i).z() += disp;

	}
}
//=========================================================================
void RBF::ResetMesh(meshPtr_Type& mesh, const meshPtr_Type& mesh_ref )
{
	for (int i = 0; i < mesh->numVertices(); i++)
	{
		mesh->point(i).x() = mesh_ref->point(i).x();

		mesh->point(i).y() = mesh_ref->point(i).y();

		mesh->point(i).z() = mesh_ref->point(i).z();
	}
}
//=========================================================================
void RBF::ExportReferencePoints(std::string FilePrefix)
{

 	std::ofstream myfile(FilePrefix+"_RBFref.csv");

 	myfile << "x coord, y coord, z coord\n";
 	for (int i = 0; i < M_NumPoints; i++)
 	{

 		myfile << std::fixed << std::setprecision(5) << M_InterpPoints[i][0] << ", ";
 		myfile << std::fixed << std::setprecision(5) << M_InterpPoints[i][1] << ", ";
 		myfile << std::fixed << std::setprecision(5) << M_InterpPoints[i][2];

 		myfile << "\n";
 	}

	myfile.close();

}
//=========================================================================
void RBF::ExportDeformedPoints(std::string FilePrefix)
{
	std::ofstream myfile(FilePrefix+"_RBFdef.csv");

	myfile << "x coord, y coord, z coord\n";
	for (int i = 0; i < M_NumPoints; i++)
	{

		myfile << std::fixed << std::setprecision(5) << M_Displacement_InterpPoints[i][0]+M_InterpPoints[i][0] << ", ";
		myfile << std::fixed << std::setprecision(5) << M_Displacement_InterpPoints[i][1]+M_InterpPoints[i][1] << ", ";
		myfile << std::fixed << std::setprecision(5) << M_Displacement_InterpPoints[i][2]+M_InterpPoints[i][2];

		myfile << "\n";
	}

	myfile.close();

}
//=========================================================================
