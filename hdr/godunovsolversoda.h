#pragma once

#include "abstractsolver.h"
struct GodunovSolverSoda: public AbstractSolver
{
    GodunovSolverSoda(Mixture mixture_, macroParam startParam_, solverParams solParam_, SystemOfEquationType type,RiemannSolverType riemannType):
        AbstractSolver(mixture_,startParam_,solParam_, type,riemannType){};

    // запускает процесс решения задачи
    void solve();

    // устанавливает некоторые граничные условия (TODO сделать более общую структуру)
    void setBorderConditions(BorderConditionSoda borderSoda_){borderSoda = borderSoda_;};

    // уникальные граничные условия
    BorderConditionSoda borderSoda;
protected:

    //заполняет начальные ячеки
    void prepareSolving() ;

    void updatePoints();


    // обновляет вектор U
    void updateU();


    double shareViscositySimple(double currentT);


    vector<macroParam>rezultAfterPStart;
};
