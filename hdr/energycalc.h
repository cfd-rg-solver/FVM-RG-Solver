#pragma once
#include "macroparam.h"

struct EnergyCalc
{
    EnergyCalc(){};
    virtual double calcEnergy(macroParam & point) = 0;
    auto getEnergyFunc(){
        return [&](macroParam & point) -> double{
            return this->calcEnergy(point);
        };
    }
    virtual double getEntalpTotal(macroParam & point){return getEntalp(point,0) + getEntalp(point,1);};
    virtual double getEntalp(macroParam & point, size_t component) = 0;
};
struct OneTempApprox : public EnergyCalc
{
    OneTempApprox(){};
    double calcEnergy(macroParam & point);
    double getEntalpTotal(macroParam & point){return getEntalp(point,0) + getEntalp(point,1);};
    double getEntalp(macroParam & point, size_t component);

private:
    double getTrRotEnegry(macroParam & point, size_t component);
    double getVibrEnergy(macroParam & point, size_t component);
    double avgVibrEnergyDiff(macroParam &point, size_t component);
    double avgVibrEnergy(macroParam &point, size_t component);
    double vibrEnergyLvl(int lvl, macroParam &point, size_t component);
    double ZvibrDiff(macroParam &point, size_t component);
    double Zvibr(macroParam &point, size_t component);
};
