#include <taurus/Core/Utilities/VITreader.hpp>

using namespace LifeV;

VITreader::VITreader() :
	M_numConfigurations (0),
	M_numParameters  (0)
{}

VITreader::~VITreader()
{}

void VITreader::Setup(const std::string FileName)
{
	std::fstream fin(FileName);
	if (!fin.is_open())
	    std::perror("Error while opening the vit file");

	std::string line;

	int check = 0;
	int index = -1;

	while(std::getline(fin, line))
	{
		if ( line[0] != '%' && check == 0 )
		{
			std::stringstream(line) >> M_numConfigurations >> M_numParameters;
			M_parameters.resize(M_numConfigurations);
			check = 1;
		}
		else if( line[0] != '%' && check == 1 )
		{
			++index;
			std::istringstream buffer(line);
			M_parameters[index].resize(M_numParameters);
			std::vector<double> myline((std::istream_iterator<double>(buffer)), std::istream_iterator<double>());
			M_parameters[index] = myline;
		}
	}
}

void VITreader::showMe( const boost::shared_ptr< Epetra_Comm >& comm ) const
{
	if ( comm->MyPID() == 0 )
	{
		std::cout << "\n\nNumber of Configurations: " << M_numConfigurations
				  << ", Number of Parameters: " << M_numParameters << std::endl;

		for (int i = 0; i < M_numConfigurations; i++)
		{
			std::cout << "Configuration " << i+1 << ": ";
			for (int j = 0; j < M_numParameters; j++)
			{
				std::cout << M_parameters[i][j] << "  ";
			}
			std::cout << "\n";
		}
	}
}
