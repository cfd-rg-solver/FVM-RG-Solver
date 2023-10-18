#include "coeffsolver.h"
#include "global.h"
#include <cmath>

double CoeffSolver::shareViscositySimple(macroParam currentPoint)
{
    double temp1 = sqrt(kB*currentPoint.temp/(M_PI*currentPoint.mixture.components[0].mass));
    double temp2 = 2.*M_PI*pow(currentPoint.mixture.components[0].sigma,2);
    double omega2 = temp1*temp2;
    return (5.*kB*currentPoint.temp) /(8.*omega2);
}
double CoeffSolver::shareViscositySimple(macroParam currentPoint, double temperature)
{
    double temp1 = sqrt(kB*temperature/(M_PI*currentPoint.mixture.components[0].mass));
    double temp2 = 2.*M_PI*pow(currentPoint.mixture.components[0].sigma,2);
    double omega2 = temp1*temp2;
    return (5.*kB*temperature) /(8.*omega2);
}
double CoeffSolver::lambda(macroParam currentPoint)
{
    double temp1 = sqrt(kB*currentPoint.temp/(M_PI*currentPoint.mixture.components[0].mass));
    double temp2 = 2.*M_PI*pow(currentPoint.mixture.components[0].sigma,2);
    double omega2 = temp1*temp2;
    return (75.*pow(kB,2)*currentPoint.temp) /(32. * currentPoint.mixture.components[0].mass * omega2 ) ;
}

double CoeffSolver::lambda(macroParam currentPoint, double temperature)
{
    double temp1 = sqrt(kB*temperature/(M_PI*currentPoint.mixture.components[0].mass));
    double temp2 = 2.*M_PI*pow(currentPoint.mixture.components[0].sigma,2);
    double omega2 = temp1*temp2;
    return (75.*pow(kB,2)*temperature) /(32. * currentPoint.mixture.components[0].mass * omega2 ) ;
}

double CoeffSolver::shareViscosityOmega(Mixture mix,double currentT)
{
    return (5.*kB*currentT) /(8.*getOmega22(mix,currentT));
}

double CoeffSolver::getOmega22(Mixture mix, double T)
{
    vector<double> f = {-0.40811, -0.05086, 0.34010, 0.70375, -0.10699, 0.00763};
    double a22= 1.5;
    double x = (log(T/mix.epsilonDevK(0)))+a22; //! Mistake in data for Ar

    double omegaLD = pow(f[0] + f[1]/pow(x,2) + f[2]/x + f[3]*x + f[4]*pow(x,2)+f[5]*pow(x,3),-1);
    double omegaS = sqrt(kB*T/(M_PI*mix.mass(0)))*2.*M_PI*pow(mix.sigma(0),2);
    return omegaLD*omegaS;
}

double CoeffSolver::bulcViscositySimple(Mixture mix,double currentT, double density, double pressure)
{
    //double Nav = 6.02214129e23;
    double n = density/mix.mass(0);
    double Crot = kB/mix.mass(0);
    double Ctr = 3.0*Crot/2.0; //! Mistake
    double Cu = Crot + Ctr;
    double F = 1+ pow(M_PI,3./2.)/2.*pow(kB*currentT/2.,-1/2.) + (pow(M_PI,2)/4. +2.)*pow(currentT/mix.epsilonDevK(0),-1) + pow(M_PI,3./2.)*pow(currentT/mix.epsilonDevK(0),-3./2.);
    double ZettaRot = ZettaInf/F;
    double tauR = (ZettaRot*M_PI*shareViscosityOmega(mix, currentT)/(4.*pressure));
    double Brot = (3.*Crot)/(2.*n*Cu*tauR);
    // return (kB*currentT/Brot)* pow(Crot/Cu),2); //! Mistake, internal specific heat for Ar is 0, and here it is obviously not...
    return 0;
}
