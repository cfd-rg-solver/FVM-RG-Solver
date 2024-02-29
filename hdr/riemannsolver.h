#pragma once

#include "global.h"
#include "systemofequation.h"
enum RiemannSolverType
{
    HLLCSolver,         // doesn't work
    HLLESolver,         // doesn't work
    HLLSimple,          // doesn't work
    ExacRiemanSolver,   // doesn't work

};
struct RiemannSolver
{
    RiemannSolver(){};
    double maxSignalVelocity = 0;
    void toMaxVelocity(double vel); // если ввести -1, то значение максимальной сокрости обнулится
    virtual void computeFlux(SystemOfEquation *system){};
    virtual void computeFlux(SystemOfEquation *system, double dt, double dh){};
    virtual void computeFlux(SystemOfEquation *system, double dh){};

    solverParams solParam;
};

struct HLLCSolver : public RiemannSolver
{
    HLLCSolver(){};
    void computeFlux(SystemOfEquation *system);
    void computeFlux(SystemOfEquation *system, double dt, double dh);
};


struct HLLESolver : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system);
};


struct HLLSimple : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system, double dt, double dh);
    void computeFlux(SystemOfEquation *system);
};

struct ExacRiemanSolver : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system, double dh);

private:
     macroParam exacRiemanSolver(macroParam left, macroParam right, double Gamma);
};
