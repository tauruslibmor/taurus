#include "Taurus_config.h"

#if defined(HAVE_LIBGP) && defined(HAVE_EIGEN)

#include <taurus/Core/Utilities/GPmodel.hpp>

using namespace std;

//=========================================================================
// Constructor //
GPmodel::GPmodel() :
				M_NumTrainingPoints  ( 0 )
{
    M_gp = new GaussianProcess (1, "CovSum ( CovSEiso, CovNoise)");
    M_hyperParams.resize(3);
    M_hyperParams << 2.0, 2.0, -1.5;
    M_gp->covf().set_loghyper(M_hyperParams);
    
    // print out initial hyperparameters
    cout << "\n Hyp0: loghyper parameter 0 = " << M_gp->covf().get_loghyper()(0);
    cout << "\n Hyp0: loghyper parameter 1 = " << M_gp->covf().get_loghyper()(1);
    cout << "\n Hyp0: loghyper parameter 2 = " << M_gp->covf().get_loghyper()(2) << "\n";

}
//=========================================================================
void GPmodel::AddTrainingPoint(const double& x, const double& y)
{
    double tmp[] = {x};
    M_gp->add_pattern(tmp, y);
    
    M_TrainingPoints_x.push_back(x);
    M_TrainingPoints_y.push_back(y);
    M_NumTrainingPoints = M_NumTrainingPoints + 1;
}
//=========================================================================
void GPmodel::SetInitHyperParameters(const double& p1, const double& p2, const double& p3)
{
    M_hyperParams.resize(3);
    M_hyperParams << p1, p2, p3;
    M_gp->covf().set_loghyper(M_hyperParams);
    
    // print out initial hyperparameters
    cout << "\n Hyp0: loghyper parameter 0 = " << M_gp->covf().get_loghyper()(0);
    cout << "\n Hyp0: loghyper parameter 1 = " << M_gp->covf().get_loghyper()(1);
    cout << "\n Hyp0: loghyper parameter 2 = " << M_gp->covf().get_loghyper()(2) << "\n";
}
//=========================================================================
void GPmodel::Generate()
{
    // Optimize hyper parameters
    int numOptIter = 150;
    bool verbosityOpt = true;
    double eps_stop = 0;
    
    M_optimizer.init(eps_stop);
    M_optimizer.maximize(M_gp, numOptIter, verbosityOpt);
    
    std::cout << "\n Optimized loghyper parameter 0 = " << M_gp->covf().get_loghyper()(0);
    std::cout << "\n Optimized loghyper parameter 1 = " << M_gp->covf().get_loghyper()(1);
    std::cout << "\n Optimized loghyper parameter 2 = " << M_gp->covf().get_loghyper()(2) << "\n";
}
//=========================================================================
void GPmodel::Save(const string& filename)
{
    M_gp->write(filename.c_str());
}
//=========================================================================
void GPmodel::Evaluate( const double& x, const double& conf_level, double& mean, double& upper_conf, double& lower_conf)
{
    double s2n = exp(2*M_gp->covf().get_loghyper()(2));

    double xx[] = {x};
    mean = M_gp->f(xx);
    upper_conf = mean + 1.96*std::sqrt(M_gp->var(xx)+s2n);
    lower_conf = mean - 1.96*std::sqrt(M_gp->var(xx)+s2n);
}
//=========================================================================
void GPmodel::Evaluate( const std::vector<double>& x, const string filename)
{
    int NumTestingPoints = x.size();
    std::ofstream myTfile;
    myTfile.open (filename);
    myTfile << "T_set = [";
    
    double f, f_up, f_down;
    for(int i = 0; i < NumTestingPoints; ++i)
    {
        Evaluate(x[i], 0.05, f, f_up, f_down);
        myTfile << x[i] << " " << f << " " << f_up << " " << f_down << "\n" << std::flush;
    }
    myTfile << "];";
    
    myTfile.close();
}
//=========================================================================
void GPmodel::Evaluate( int NumTestingPoints, const string filename)
{
    // extreme values
    double min_x = *std::min_element(  M_TrainingPoints_x.begin(), M_TrainingPoints_x.end() );
    double max_x = *std::max_element(  M_TrainingPoints_x.begin(), M_TrainingPoints_x.end() );
    
    min_x = min_x*0.95;
    max_x = max_x*1.05;
    
    double h = (max_x - min_x) / NumTestingPoints;
    
    std::vector<double> test_points;
    
    for(int i = 0; i < NumTestingPoints+1; ++i)
    {
        test_points.push_back(min_x+i*h);
    }
    Evaluate( test_points, filename);
}
//=========================================================================
void GPmodel::WriteTraining(const string& filename)
{
    std::ofstream myTfile;
    myTfile.open (filename);
    myTfile << M_NumTrainingPoints << std::flush;;
    
    double f, f_up, f_down;
    for(int i = 0; i < M_NumTrainingPoints; ++i)
    {
        myTfile << "\n" << M_TrainingPoints_x[i] << " " << M_TrainingPoints_y[i] << std::flush;
    }
    
    myTfile.close();
}
//=========================================================================
void GPmodel::LoadTraining(const string& filename)
{
    // read from file training points
    std::fstream myfile(filename, std::ios_base::in);
    
    myfile >> M_NumTrainingPoints;
    std::cout << "\n NumTrainingPoints = " << M_NumTrainingPoints << std::endl;
    
    M_TrainingPoints_x.resize(M_NumTrainingPoints);
    M_TrainingPoints_y.resize(M_NumTrainingPoints);
    
    for (int i = 0; i < M_NumTrainingPoints; i++)
    {
        myfile >> M_TrainingPoints_x[i];
        myfile >> M_TrainingPoints_y[i];
        double x[] = {M_TrainingPoints_x[i]};
        M_gp->add_pattern(x, M_TrainingPoints_y[i]);
    }
    
    myfile.close();
}
//=========================================================================
#endif
    

