#pragma once
#include "macroparam.h"
#include "energycalc.h"
#include <functional>

struct NonLinearEqSolver
{
	virtual double solveEq(EnergyCalc* energy, macroParam point, double rightPart) = 0;
};
struct Newton : public NonLinearEqSolver
{
	double solveEq(EnergyCalc* energy, macroParam point, double rightPart);
};

struct IntegralCalculator
{
	virtual double integrate(std::vector<double> x, std::vector<double> y) = 0;
	virtual double integrate(std::vector<double> y, double dh) = 0;
};

struct IntegralCalculatorTrap : public IntegralCalculator
{
	double integrate(std::vector<double> x, std::vector<double> y);
	double integrate(std::vector<double> y, double dh);
};

struct DerrivativeCalculator
{
	static double derrivativeCenter(std::vector<double> y, double dh, int idx);
	static double derrivativeLeft(std::vector<double> y, double dh, int idx);
	static double derrivativeRight(std::vector<double> y, double dh, int idx);
};