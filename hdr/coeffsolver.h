#pragma once
#include "global.h"

struct CoeffSolver
{
public:
    virtual double shareViscositySimple(macroParam currentPoint) = 0;
    virtual double lambda(macroParam currentPoint) = 0;
    virtual double shareViscosityOmega(Mixture mix,double currentT) = 0;
    virtual double getOmega22(Mixture mix,double T) = 0;
    virtual double bulcViscositySimple(macroParam currentPoint) = 0;
    virtual double binaryDiffusion(macroParam currentPoint){return 0;};
    virtual double omega11(macroParam currentPoint){return 0;};

};
struct CoeffSolver1Comp1Temp : public CoeffSolver
{
    double shareViscositySimple(macroParam currentPoint);
    double lambda(macroParam currentPoint);
    double shareViscosityOmega(Mixture mix,double currentT);
    double getOmega22(Mixture mix,double T);
    double bulcViscositySimple(macroParam currentPoint);
};
struct CoeffSolver2Comp1Temp : public CoeffSolver1Comp1Temp
{
    double shareViscositySimple(macroParam currentPoint);
    double lambda(macroParam currentPoint);
    double shareViscosityOmega(Mixture mix,double currentT);
    double getOmega22(Mixture mix,double T);
    double bulcViscositySimple(macroParam currentPoint);
    double binaryDiffusion(macroParam currentPoint);
    double omega11(macroParam currentPoint);
private:
    double shareViscosity(macroParam currentPoint, size_t component);
    double lambda(macroParam currentPoint, size_t component);
    double phi(macroParam currentPoint, size_t component1, size_t component2);
};
