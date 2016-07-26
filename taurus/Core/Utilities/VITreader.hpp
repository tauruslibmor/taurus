#ifndef VITreader_H
#define VITreader_H

#include "Epetra_MpiComm.h"
#include <string>
#include <fstream>

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

class VITreader
{
	typedef std::vector<std::vector<double>> parameterContainer_Type;

public:

	VITreader();

	~VITreader();

	void Setup(const std::string FileName);

	void showMe( const boost::shared_ptr< Epetra_Comm >& comm ) const ;

	parameterContainer_Type getParameters() const { return M_parameters; };

private:

    // members
	unsigned int M_numConfigurations;
	unsigned int M_numParameters;
	parameterContainer_Type M_parameters;
};

} // namespace LifeV

#endif // VITreader_H
