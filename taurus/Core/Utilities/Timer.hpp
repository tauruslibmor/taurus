#ifndef TIMER_H
#define TIMER_H

#include "Epetra_MpiComm.h"
#include <Epetra_Time.h>
#include <lifev/core/LifeV.hpp>

#include <string>
#include <fstream>


class Timer
{
	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;
	typedef boost::shared_ptr< Epetra_Time > timePtr_Type;

public:

	Timer(commPtr_Type& comm, double scale=1, std::string time_unit="s");
	~Timer();

	void setSerial();
	void StartTimer();
	double StopTimer(const std::string name);
	double StopTimer();
	void StartGlobalTimer();
	double StopGlobalTimer(const std::string name);
	double StopGlobalTimer();
    void verbose(bool verbosity=true);
    void setPrecision(const int precision);

private:

    // members
	timePtr_Type M_GlobalTimer;
	timePtr_Type M_Timer;
	commPtr_Type M_comm;
	double M_scale;
	std::string M_time_unit;
	bool M_verbose;
	bool M_isSerial;
	int M_precision;

	//methods

	double maxTime(double local_time);
	double minTime(double local_time);
	double averageTime(double local_time);
	double ConvertTime(double time);
};

#endif // TIMER_H
