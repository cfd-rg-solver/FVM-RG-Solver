#pragma once
#include "global.h"

struct CoeffSolver
{
public:

    double shareViscositySimple(macroParam currentPoint);
    double shareViscositySimple(macroParam currentPoint, double temperature);
    double lambda(macroParam currentPoint);
    double lambda(macroParam currentPoint, double temperature);
    double shareViscosityOmega(Mixture mix,double currentT);
    double getOmega22(Mixture mix,double T);
    double bulcViscositySimple(Mixture mix,double currentT, double density, double pressure);

};
