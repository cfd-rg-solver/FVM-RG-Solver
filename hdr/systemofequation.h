#pragma once
#include "global.h"
#include "coeffsolver.h"
#include "bordercondition.h"
#include "energycalc.h"
#include "numeric.h"
enum SystemOfEquationType
{
    couette2,
    couette2Alt,
    couette2AltBinary,
    soda
};
struct SystemOfEquation
{
    SystemOfEquation(){};
    void setNumberOfCells( size_t cells);
    void setMixture(Mixture mixture_);
    void setBorderCondition(BorderCondition* border_);
    void setCoeffSolver(CoeffSolver* coeffSolver_);
    void setSolverParams(solverParams solParam_);
    void setEqSolver(NonLinearEqSolver* eqSolver_);

    virtual double getPressure(size_t i) = 0;
    virtual double getDensity(size_t i) = 0;
    virtual double getDensity(size_t i, size_t component);
    virtual double getVelocity(size_t i) = 0;
    virtual double getVelocityTau(size_t i) = 0;
    virtual double getVelocityNormal(size_t i) = 0;
    virtual double getEnergy(size_t i) = 0;
    virtual double getTemp(size_t i) = 0;

    double getMaxVelocity();

    virtual void prepareIndex() = 0;
    virtual void prepareVectorSizes();
    virtual void prepareSolving(vector<macroParam> & points) = 0;

    virtual void updateU(double dh, double dt) = 0;
    virtual void updateBorderU(vector<macroParam> & points) = 0;
    virtual void computeF(vector<macroParam> & points, double dh) = 0;
    virtual void computeFv(vector<macroParam> & points, double dh){};

    size_t numberOfCells;
    size_t numberOfComponents;

    // индексы соответсвующих компонент в векторах U,F,R
    size_t v_tau;
    size_t v_normal;
    size_t energy;

    // сколько всего уравнений (не индекс)
    size_t systemOrder;

    Mixture mixture;
    EnergyCalc* energyCalculator = new OneTempApprox();
    NonLinearEqSolver*  eqSolver;
    CoeffSolver* coeffSolver;
    BorderCondition* border;
    solverParams solParam;
    SystemOfEquationType systemType;

    vector<Matrix> U, R, F, Flux;
};


struct Couette2 : public SystemOfEquation
{
    Couette2(){systemType = SystemOfEquationType::couette2;};
    void prepareSolving(vector<macroParam> & points);
    void prepareIndex();

    double getPressure(size_t i);
    double getDensity(size_t i);
    double getVelocity(size_t i);
    double getVelocityTau(size_t i);
    double getVelocityNormal(size_t i);
    double getEnergy(size_t i);
    double getTemp(size_t i);

    void updateU(double dh, double dt);
    void updateBorderU(vector<macroParam> & points);
    void computeF(vector<macroParam> & points, double dh);

};

struct Couette2Alt : public Couette2
{
    Couette2Alt(){systemType = SystemOfEquationType::couette2Alt;};

    void prepareVectorSizes();

    void updateU(double dh, double dt);
    void computeF(vector<macroParam> & points, double dh);
    void computeFv(vector<macroParam> & points, double dh);

    vector<Matrix> Fv;
};

struct Couette2AltBinary : public Couette2Alt
{
    Couette2AltBinary(){systemType = SystemOfEquationType::couette2AltBinary;};

    void prepareSolving(vector<macroParam> & points);

    double getPressure(size_t i);
    double getDensity(size_t i);
    double getDensity(size_t i, size_t component);
    double getTemp(size_t i);

    void updateU(double dh, double dt);
    void updateBorderU(vector<macroParam> & points);
    void computeF(vector<macroParam> & points, double dh);
    void computeFv(vector<macroParam> & points, double dh);
private:
    std::vector<double> temperature;
    void calcAndRemeberTemp();
};

struct Soda : public SystemOfEquation
{
    Soda(){systemType = SystemOfEquationType::soda;};
    void prepareSolving(vector<macroParam> & points);
    void prepareIndex();

    double getPressure(size_t i);
    double getDensity(size_t i);
    double getVelocity(size_t i);
    double getVelocityTau(size_t i);
    double getVelocityNormal(size_t i);
    double getSoundSpeed(size_t i);
    double getEnergy(size_t i);
    double getTemp(size_t i) { return 0; };

    double getMaxVelocity();
    void updateU(double dh, double dt);
    void updateBorderU(vector<macroParam> & points){};
    void computeF(vector<macroParam> & points, double dh);


    double gamma = 1.4;
};
