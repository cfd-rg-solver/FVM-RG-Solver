#include "coeffsolver.h"
#include "global.h"
#include "energycalc.h"
#include <cmath>

double CoeffSolver1Comp1Temp::shareViscositySimple(macroParam currentPoint)
{
    double temp1 = sqrt(kB*currentPoint.temp/(M_PI*currentPoint.mixture.components[0].mass));
    double temp2 = 2.*M_PI*pow(currentPoint.mixture.components[0].sigma,2);
    double omega2 = temp1*temp2;
    return (5.*kB*currentPoint.temp) /(8.*omega2);
}

double CoeffSolver1Comp1Temp::lambda(macroParam currentPoint)
{
    double temp1 = sqrt(kB*currentPoint.temp/(M_PI*currentPoint.mixture.components[0].mass));
    double temp2 = 2.*M_PI*pow(currentPoint.mixture.components[0].sigma,2);
    double omega2 = temp1*temp2;
    return (75.*pow(kB,2)*currentPoint.temp) /(32. * currentPoint.mixture.components[0].mass * omega2 ) ;
}

double CoeffSolver1Comp1Temp::lambdaMultiAtom(macroParam currentPoint)
{
    size_t component = 0;  // consider 1 component gas

    double lambda_tr = (75. * pow(kB, 2) * currentPoint.temp) / (32. * currentPoint.mixture.mass(component) * getOmega22(currentPoint.mixture, currentPoint.temp));

    double Crot = 3. / 2. * kB / currentPoint.mixture.mass(component); // because methane has 3 degrees of freedom
    // OneTempApproxMultiModes OneTempApprox;
    double Cvibr = OneTempApprox.getCvibr(currentPoint, component);
    double Cint = Crot + Cvibr;
    double Omega11 = getOmega11(currentPoint.mixture, currentPoint.temp);
    double lambda_int = (3. * kB * currentPoint.temp) / (8. * Omega11) * Cint;

    return lambda_tr + lambda_int;
}

double CoeffSolver1Comp1Temp::shareViscosityOmega(Mixture mix,double currentT)
{
    return (5.*kB*currentT) /(8.*getOmega22(mix,currentT));
}

double CoeffSolver1Comp1Temp::getOmega11(Mixture mix, double T)
{
    size_t component = 0; // consider 1 component gas
    double omegaRS = sqrt(kB * T / (M_PI * mix.mass(component))) * M_PI * pow(mix.sigma(component), 2);
    vector<double> f = {-0.16845, -0.02258, 0.19779, 0.64373, -0.092679, 0.00711};
    double a11= 1.4;
    double x = (log(T/mix.epsilonDevK(component)))+a11;
    double omegaLD = pow(f[0] + f[1]/pow(x,2) + f[2]/x + f[3]*x + f[4]*pow(x,2)+f[5]*pow(x,3),-1);
    return omegaRS * omegaLD;
}

double CoeffSolver1Comp1Temp::getOmega22(Mixture mix, double T)
{
    vector<double> f = {-0.40811, -0.05086, 0.34010, 0.70375, -0.10699, 0.00763};
    double a22= 1.5;
    double x = (log(T/mix.epsilonDevK(0)))+a22; //! Mistake in data for Ar

    double omegaLD = pow(f[0] + f[1]/pow(x,2) + f[2]/x + f[3]*x + f[4]*pow(x,2)+f[5]*pow(x,3),-1);
    double omegaS = sqrt(kB*T/(M_PI*mix.mass(0)))*2.*M_PI*pow(mix.sigma(0),2);
    return omegaLD*omegaS;
}

double CoeffSolver1Comp1Temp::bulcViscositySimple(macroParam currentPoint)
{
    //double Nav = 6.02214129e23;
    double density = currentPoint.density;
    double m = currentPoint.mixture.mass(0);
    double n = density/m;
    double T = currentPoint.temp;
    double p = currentPoint.pressure;
    double epsilon = currentPoint.mixture.epsilonDevK(0);

    double Crot = kB/m;
    double Ctr = 3.0*Crot/2.0;
    double Cu = Crot + Ctr;
    double F = 1+ pow(M_PI,3./2.)/2.*pow(kB*T/2.,-1/2.) + (pow(M_PI,2)/4. +2.)*pow(T/epsilon,-1) + pow(M_PI,3./2.)*pow(T/epsilon,-3./2.);
    double ZettaRot = ZettaInf/F;
    double tauR = (ZettaRot*M_PI*shareViscosityOmega(currentPoint.mixture, T)/(4.*p));
    double Brot = (3.*Crot)/(2.*n*Cu*tauR);
    // return (kB*currentT/Brot)* pow(Crot/Cu),2); //! Mistake, internal specific heat for Ar is 0, and here it is obviously not...
    return 0;
}

double CoeffSolver1Comp1Temp::bulkViscosityMultiAtom(macroParam point)
{
    /* function to calculate bulk viscosity zeta = kB*T/beta_int * (c_int/c_V)^2 */

    size_t component = 0; // consider 1 component gas

    double Ctr = 3. / 2. * kB / point.mixture.mass(component);
    double Crot = 3. / 2. * kB / point.mixture.mass(component); // because methane has 3 degrees of freedom
    // OneTempApproxMultiModes OneTempApprox;
    double Cvibr = OneTempApprox.getCvibr(point, component);
    double Cv = Ctr + Crot + Cvibr;
    double Cint = Cv - Ctr;

    double Zinf = point.mixture.Zinf(component);
    double eps = point.mixture.epsilonDevK(component) * kB;
    double F = 1 + pow(M_PI, 3. / 2.) / 2. * pow(kB * point.temp / eps, -1. / 2.) + (pow(M_PI, 2) / 4. + 2.) * pow(kB * point.temp / eps, -1) + pow(M_PI, 3. / 2.) * pow(kB * point.temp / eps, -3. / 2.);
    double Zrot = Zinf / F;
    double eta = shareViscosityOmega(point.mixture, point.temp);

    double tau_rot = M_PI * eta * Zrot / (4. * point.pressure); // Parker model

    double tau_vibr = exp(-5.4 + 40. * pow(point.temp, -1. / 3.)) * 10e-6; // Vibrational Relaxation of METHANE L.Willard Richards and David H.Sigafoos

    double beta_int = point.pressure * UniversalGasConstant * pow(Cint / Cv, 2) / (Crot / tau_rot + Cvibr / tau_vibr);

    double zeta = kB * point.temp * pow(Cint / Cv, 2) / (beta_int);
    return zeta;
}

double CoeffSolver2Comp1Temp::shareViscositySimple(macroParam currentPoint)
{
    double m1 = currentPoint.mixture.mass(0);
    double M1 = currentPoint.mixture.molarMass(0);
    double y1 =  currentPoint.fractionArray[0];
    double etta1 = shareViscosity(currentPoint,0);

    double m2 = currentPoint.mixture.mass(1);
    double M2 = currentPoint.mixture.molarMass(1);
    double y2 =  currentPoint.fractionArray[1];
    double etta2 = shareViscosity(currentPoint,1);

    double phi12 = phi(currentPoint,0,1);
    double phi21 = phi(currentPoint,1,0);

    double etta = etta1 * pow(1 + (m1 * y2) / (m2 * y1) * phi12, -1) + etta2 * pow(1 + (m2 * y1) / (m1 * y2) * phi21, -1);
    return etta;
}

double CoeffSolver2Comp1Temp::lambda(macroParam currentPoint)
{
    double m1 = currentPoint.mixture.mass(0);
    double M1 = currentPoint.mixture.molarMass(0);
    double y1 =  currentPoint.fractionArray[0];
    double lambda1 = lambda(currentPoint,0);

    double m2 = currentPoint.mixture.mass(1);
    double M2 = currentPoint.mixture.molarMass(1);
    double y2 =  currentPoint.fractionArray[1];
    double lambda2 = lambda(currentPoint,1);

    double G12 = 1.065 * phi(currentPoint,0,1);
    double G21 = 1.065 * phi(currentPoint,1,0);

    double lambda = lambda1 * pow(1 + (m2 * y1) / (m1 * y2) * G12, -1)  +  lambda2 * pow(1 + (m1 * y2) / (m2 * y1) * G21, -1);
    return lambda;
}

double CoeffSolver2Comp1Temp::bulcViscositySimple(macroParam currentPoint)
{
    return 0;
}

double CoeffSolver2Comp1Temp::effDiffusion(macroParam currentPoint, size_t component)
{
    if(currentPoint.mixture.NumberOfComponents == 1)
        return 0;


    double xComp = Xc(currentPoint,component);
    double denominator = 0;
    for(int i = 0 ; i < currentPoint.mixture.NumberOfComponents; i++)
    {
        denominator += Xc(currentPoint,i)/binaryDiffusion(currentPoint,component,i);
    }
    double Dcomp = (1 - xComp) / denominator;
    return Dcomp;
}

double CoeffSolver2Comp1Temp::binaryDiffusion(macroParam currentPoint, size_t comp1, size_t comp2)
{
    double M = 0;
    if(comp1 != comp2)
        M = 1 / ( currentPoint.fractionArray[comp1] / currentPoint.mixture.molarMass(comp1) +
                 currentPoint.fractionArray[comp2] / currentPoint.mixture.molarMass(comp2));
    else
        M = currentPoint.mixture.molarMass(comp1);
    double R = UniversalGasConstant;
    double m1 = currentPoint.mixture.mass(comp1);
    double m2 = currentPoint.mixture.mass(comp2);
    double D21 = (3 * pow(kB,2) * M * currentPoint.temp * (m1 + m2)) / (16 * currentPoint.density * R * m1 * m2 ) / getOmega11(currentPoint,comp1,comp2);
    return D21;
}

double CoeffSolver2Comp1Temp::Xc(macroParam currentPoint, size_t component)
{
    double rhoComp = currentPoint.densityArray[component];
    double Mcomp = currentPoint.mixture.molarMass(component);

    double n = 0;
    for(int i = 0 ; i < currentPoint.mixture.NumberOfComponents; i++)
    {
        double tmp = currentPoint.densityArray[i];
        for(int j = 0 ; j < currentPoint.mixture.NumberOfComponents; j++)
        {
            if(j == i)
                continue;
            tmp *=  currentPoint.mixture.molarMass(j);
        }
        n += tmp;
    }
    double xComp = (rhoComp * Mcomp) / n;
    return xComp;
}

double CoeffSolver2Comp1Temp::shareViscosity(macroParam currentPoint, size_t component)
{
    double m = currentPoint.mixture.components[component].mass;
    double sigma = currentPoint.mixture.components[component].sigma;
    double temp1 = sqrt(kB*currentPoint.temp/(M_PI*m));
    double temp2 = 2.*M_PI*pow(sigma,2);
    double omega2 = temp1*temp2;
    return (5.*kB*currentPoint.temp) /(8.*omega2);
}

double CoeffSolver2Comp1Temp::lambda(macroParam currentPoint, size_t component)
{
    double m = currentPoint.mixture.components[component].mass;
    double sigma = currentPoint.mixture.components[component].sigma;

    double temp1 = sqrt(kB*currentPoint.temp/(M_PI*m));
    double temp2 = 2.*M_PI*pow(sigma,2);
    double omega2 = temp1*temp2;
    return (75.*pow(kB,2)*currentPoint.temp) /(32. * m * omega2 ) ;
}

double CoeffSolver2Comp1Temp::phi(macroParam currentPoint, size_t component1, size_t component2)
{
    double M1 = currentPoint.mixture.molarMass(component1);
    double etta1 = shareViscosity(currentPoint,component1);

    double M2 = currentPoint.mixture.molarMass(component2);
    double etta2 = shareViscosity(currentPoint,component2);

    double phi12 = (1 + sqrt(etta1/etta2)*pow(M2/M1,0.25)) / (2 * sqrt(2) * sqrt( 1 +  M1/M2 ));
    return phi12;
}
double CoeffSolver2Comp1Temp::getOmega11(macroParam currentPoint, size_t comp1, size_t comp2)
{
    double tmp1 = sqrt(kB * currentPoint.temp / (M_PI * (currentPoint.mixture.mass(comp1) + currentPoint.mixture.mass(comp2))));
    double sigma = ( currentPoint.mixture.sigma(comp1) + currentPoint.mixture.sigma(comp2)) / 2. ;
    double omega11 = tmp1 * M_PI * pow(sigma,2);
    return omega11;
}
