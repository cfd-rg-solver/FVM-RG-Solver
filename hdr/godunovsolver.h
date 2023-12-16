#pragma once

#include "abstractsolver.h"
struct GodunovSolver: public AbstractSolver
{
    GodunovSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_, SystemOfEquationType type,RiemannSolverType riemannType):
        AbstractSolver(mixture_,startParam_,solParam_, type,riemannType){};
    GodunovSolver(Mixture mixture_, vector<macroParam> startParam_, solverParams solParam_, SystemOfEquationType type,RiemannSolverType riemannType):
        AbstractSolver(mixture_,startParam_,solParam_, type,riemannType){};


    // запускает процесс решения задачи
    void solve();
protected:

    // Расчет релаксационных членов
    void computeR();

    // обновляет вектор U
    void updateU();

    vector<macroParam>rezultAfterPStart;
};
