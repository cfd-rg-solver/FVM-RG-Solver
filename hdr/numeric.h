#pragma once
#include "macroparam.h"
#include "energycalc.h"
#include <functional>

struct NonLinearEqSolver
{
    virtual double solveEq(EnergyCalc *energy, macroParam point, double rightPart) = 0;
};
struct Newton : public NonLinearEqSolver
{
    double solveEq(EnergyCalc *energy, macroParam point, double rightPart);
};
