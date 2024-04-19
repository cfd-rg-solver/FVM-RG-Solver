#pragma once
#include "mixture.h"

struct macroParam
{
    macroParam(){mixture = Mixture();}
    macroParam(Mixture mix):mixture(mix){fractionArray.resize(mix.NumberOfComponents); densityArray.resize(mix.NumberOfComponents);}
    Mixture mixture;
    vector<double> densityArray;
    vector<double> fractionArray; // в сумме должны давать 1
    double density      = 0;
    double pressure     = 0;
    double velocity_tau     = 0;
    double velocity_normal  = 0;
    double velocity = 0;
    double temp         = 0;
    double temp_vibr    = 0;
    double tempIntr     = 0;
    double soundSpeed   = 0;
    double gamma; // показатель адиабаты
    string gas         ="Ar";
};
