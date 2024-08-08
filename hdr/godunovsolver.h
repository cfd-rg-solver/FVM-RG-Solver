#pragma once
#include "lawchecker.h"
#include "abstractsolver.h"

struct GodunovSolver : public AbstractSolver
{
	GodunovSolver(Mixture mixture_, solverParams solParam_, SystemOfEquationType type, RiemannSolverType riemannType) :
		AbstractSolver(mixture_, solParam_, type, riemannType) {};


	// запускает процесс решения задачи
	void solve();

	void setOutputDirectory(const std::string& outputDir);
protected:

	// Расчет релаксационных членов
	void computeR();

	// обновляет вектор U
	void updateU();

	void printViscousFluxComponents(SystemOfEquation* system, int N);
		//    vector<macroParam>rezultAfterPStart;

	std::string outputDirectory;
};
