#pragma once
#include "global.h"
#include "macroparam.h"
#include "energycalc.h"
struct CoeffSolver
{
public:
	virtual double shareViscositySimple(macroParam currentPoint) = 0;
	virtual double lambda(macroParam currentPoint) = 0;
	virtual double lambdaMultiAtom(macroParam currentPoint) = 0;
	virtual double shareViscosityOmega(Mixture mix, double currentT) = 0;
	virtual double getOmega22(Mixture mix, double T) = 0;
	virtual double getOmega11(Mixture mix, double T) = 0;
	virtual double bulcViscositySimple(macroParam currentPoint) = 0;
	virtual double bulkViscosityMultiAtom(macroParam currentPoint) = 0;
	virtual double effDiffusion(macroParam currentPoint, size_t component) { return 0; };
	virtual double binaryDiffusion(macroParam currentPoint, size_t comp1, size_t comp2) { return 0; };
};
struct CoeffSolver1Comp1Temp : public CoeffSolver
{
	double shareViscositySimple(macroParam currentPoint);
	double lambda(macroParam currentPoint);
	double lambdaMultiAtom(macroParam currentPoint);
	double shareViscosityOmega(Mixture mix, double currentT);
	double getOmega11(Mixture mix, double T);
	double getOmega22(Mixture mix, double T);
	double bulcViscositySimple(macroParam currentPoint);
	double bulkViscosityMultiAtom(macroParam currentPoint);

private:
	OneTempApproxMultiModes OneTempApprox;
};
struct CoeffSolver2Comp1Temp : public CoeffSolver1Comp1Temp
{
	double shareViscositySimple(macroParam currentPoint);
	double lambda(macroParam currentPoint);
	double bulcViscositySimple(macroParam currentPoint);
	double effDiffusion(macroParam currentPoint, size_t component);
	double binaryDiffusion(macroParam currentPoint, size_t comp1, size_t comp2);
private:
	double Xc(macroParam currentPoint, size_t component);
	double shareViscosity(macroParam currentPoint, size_t component);
	double lambda(macroParam currentPoint, size_t component);
	double phi(macroParam currentPoint, size_t component1, size_t component2);
	double getOmega11(macroParam currentPoint, size_t comp1, size_t comp2);
};
