#include "energycalc.h"

double OneTempApprox::calcEnergy(macroParam &point)
{
    double UTrRot = getTrRotEnegry(point, 0) + getTrRotEnegry(point, 1);
    double UVibr =  getVibrEnergy(point, 0) + getVibrEnergy(point, 1);
    double E = point.density * (UTrRot + UVibr) + 0.5*pow(point.velocity,2)*point.density;
    return E;
}
double OneTempApprox::getEntalp(macroParam &point, size_t component)
{
    double Tr = 5/2. * kB * point.temp * point.fractionArray[component] / point.mixture.components[component].mass;
    double Rot = 0;
    if(point.mixture.components[component].numberAtoms == 2)
        Rot = kB * point.temp * point.fractionArray[component] / point.mixture.components[component].mass;
    double res = Tr + Rot + getVibrEnergy(point, component);
    return res;
}

double OneTempApprox::getTrRotEnegry(macroParam &point, size_t component)
{
    int i = point.mixture.components[component].numberAtoms;
    double U = (i * 2 + 1)/2. * kB * point.temp * point.fractionArray[component] / point.mixture.components[component].mass;
    return U;
}

double OneTempApprox::getVibrEnergy(macroParam &point, size_t component)
{
    if(point.mixture.components[component].numberAtoms == 1)
        return 0;
    double res = avgVibrEnergy(point , component) * point.fractionArray[component] / point.mixture.mass(component);
    return res;
}

// double OneTempApprox::getVibrEnergyMultiAtom(macroParam& point, size_t component) {
//    // todo: what changes?
//
// }

double OneTempApprox::avgVibrEnergyDiff(macroParam &point, size_t component)
{
    if(point.mixture.components[component].numberAtoms == 1)
        return 0;

    double Z = Zvibr(point, component);
    double Zdiff = ZvibrDiff(point, component);

    double tmp = 0;
    double tmpDiff = 0;
    double sum = 0;
    for(size_t i = 0; i < point.mixture.components[component].numberVibrLvl; i++)
    {
        double eps_ic = vibrEnergyLvl(i, point, component);
        tmp = eps_ic * exp (-eps_ic / ( kB * point.temp));
        tmpDiff = pow(eps_ic,2) * exp (-eps_ic / ( kB * point.temp)) / ( kB * pow(point.temp,2));
        sum += tmpDiff * Z - tmp * Zdiff;
    }

    double res = sum / (pow( Z , 2));
    return res;
}

double OneTempApprox::avgVibrEnergy(macroParam &point, size_t component)
{
    if(point.mixture.components[component].numberAtoms == 1)
        return 0;

    double sum = 0;
    double Z = Zvibr(point, component);
    for(size_t i = 0; i < point.mixture.components[component].numberVibrLvl; i++)
    {
        sum += vibrEnergyLvl(i, point, component) / Z * exp (-vibrEnergyLvl(i, point, component) / ( kB * point.temp));
    }
    return sum;
}

double OneTempApprox::vibrEnergyLvl(int lvl, macroParam &point, size_t component)
{
    return hc * point.mixture.components[component].omega_e * (lvl + 1./2.);
}

double OneTempApprox::ZvibrDiff(macroParam &point, size_t component)
{
    if(point.mixture.components[component].numberAtoms == 1)
        return 0;

    double sum = 0;
    double s_i = 1; // только для O2
    for(size_t i = 0; i < point.mixture.components[component].numberVibrLvl; i++)
    {
        double eps_ic = vibrEnergyLvl(i, point, component);
        sum += s_i * eps_ic * exp (-eps_ic / ( kB * point.temp)) / ( kB * pow(point.temp , 2.));
    }
    return sum;
}

double OneTempApprox::Zvibr(macroParam &point, size_t component)
{
    if(point.mixture.components[component].numberAtoms == 1)
        return 0;

    double sum = 0;
    double s_i = 1; // только для O2
    for(size_t i = 0; i < point.mixture.components[component].numberVibrLvl; i++)
    {
        sum += exp (-vibrEnergyLvl(i, point, component) / ( kB * point.temp));
    }
    return sum;
}
