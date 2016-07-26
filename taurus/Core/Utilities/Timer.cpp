#include <taurus/Core/Utilities/Timer.hpp>

using namespace std;

//=========================================================================
Timer::Timer(commPtr_Type& comm, double scale, std::string time_unit) :
	M_scale		   (scale),
	M_time_unit    (time_unit),
	M_comm		   (comm),
	M_verbose      (false),
	M_isSerial	   (false),
	M_precision    (3)
{
	if ( M_comm->MyPID() == 0 )
	{
		M_verbose = true;
	}

	M_GlobalTimer.reset (new Epetra_Time(*M_comm));
	M_Timer.reset (new Epetra_Time(*M_comm));
}
//=========================================================================
Timer::~Timer()
{}
//=========================================================================
void Timer::setSerial()
{
	M_isSerial = true;
	M_verbose  = true;
}
//=========================================================================
void Timer::verbose(bool verbosity)
{
    M_verbose = verbosity;
}
//=========================================================================
void Timer::setPrecision(const int precision)
{
	M_precision = precision;
}
//=========================================================================
double Timer::ConvertTime(double time)
{
	if (M_isSerial)
		return (time * M_scale);
	else
		return (maxTime(time) * M_scale);
}
//=========================================================================
double Timer::maxTime(double time)
{
	double maxtime;
	M_comm->MaxAll(&time, &maxtime, 1);
	return maxtime;
}
//=========================================================================
double Timer::minTime(double time)
{
	double mintime;
	M_comm->MinAll(&time, &mintime, 1);
	return mintime;
}
//=========================================================================
double Timer::averageTime(double time)
{
	double sumtime;
	M_comm->SumAll(&time, &sumtime, 1);
	return ( sumtime / M_comm->NumProc() );
}
//=========================================================================
void Timer::StartGlobalTimer()
{
	M_GlobalTimer->ResetStartTime();
}
//=========================================================================
void Timer::StartTimer()
{
	M_Timer->ResetStartTime();
}
//=========================================================================
double Timer::StopTimer(const std::string name)
{
	double deltaT = ConvertTime( M_Timer->ElapsedTime() );
	if (M_verbose)
	{
		std::streamsize ss = cout.precision();
		cout << "\n\n** Elapsed time for " << name << ": " << std::setprecision(M_precision) << std::scientific << deltaT << " " << M_time_unit << " ** \n\n";
		cout.unsetf(ios_base::fixed);
		cout.precision (ss);
		cout << std::fixed;
	}
	return deltaT;
}
//=========================================================================
double Timer::StopTimer()
{
	double deltaT = ConvertTime( M_Timer->ElapsedTime() );
	return deltaT;
}
//=========================================================================
double Timer::StopGlobalTimer(const std::string name)
{
	double deltaT = ConvertTime ( M_GlobalTimer->ElapsedTime() );
	if (M_verbose)
	{
		std::streamsize ss = cout.precision();
		cout << "\n\n** Elapsed time for " << name << ": " << std::setprecision(M_precision) << std::scientific << deltaT << " " << M_time_unit << " ** \n\n";
		cout.unsetf(ios_base::fixed);
		cout.precision(ss);
		cout << std::fixed;
	}
	return deltaT;
}
//=========================================================================
double Timer::StopGlobalTimer()
{
	double deltaT = ConvertTime ( M_GlobalTimer->ElapsedTime() );
	return deltaT;
}
//=========================================================================
