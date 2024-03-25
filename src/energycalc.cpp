#include "energycalc.h"
#include <iostream>

double OneTempApprox::calcEnergy(macroParam &point)
{
    double UTrRot = getTrRotEnegry(point, 0) + getTrRotEnegry(point, 1); // ! so here is binary mixture is considered
    double UVibr =  getVibrEnergy(point, 0) + getVibrEnergy(point, 1);
    double E = point.density * (UTrRot + UVibr) + 0.5*pow(point.velocity,2)*point.density; // ! so here is rho*E
    return E;
}
double OneTempApprox::getEntalp(macroParam &point, size_t component)
{
    double Tr = 5/2. * kB * point.temp * point.fractionArray[component] / point.mixture.components[component].mass; // ! why mass fraction is here?
    double Rot = 0;
    if(point.mixture.components[component].numberAtoms > 1)
        Rot = kB * point.temp * point.fractionArray[component] / point.mixture.components[component].mass;
    double res = Tr + Rot + getVibrEnergy(point, component);
    return res;
}

double OneTempApprox::getTrRotEnegry(macroParam &point, size_t component)
{
    int i = point.mixture.components[component].numberAtoms;
    double U = (i * 2 + 1)/2. * kB * point.temp * point.fractionArray[component] / point.mixture.components[component].mass; // only two-atom molecules are considered
    return U;
}

double OneTempApprox::getVibrEnergy(macroParam &point, size_t component)
{
    if(point.mixture.components[component].numberAtoms == 1)
        return 0;
    double res = avgVibrEnergy(point , component) * point.fractionArray[component] / point.mixture.mass(component); // ! here the fraction is also included
    return res;
}

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
    double s_i = 1; // ! for binary molecules
    for(size_t i = 0; i < point.mixture.components[component].numberVibrLvl; i++)
    {
        sum += s_i * vibrEnergyLvl(i, point, component) / Z * exp (-vibrEnergyLvl(i, point, component) / ( kB * point.temp));
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
    double s_i = 1; // ! for binary molecules
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
    double s_i = 1; // ! for binary molecules
    for(size_t i = 0; i < point.mixture.components[component].numberVibrLvl; i++)
    {
        sum += s_i * exp (-vibrEnergyLvl(i, point, component) / ( kB * point.temp));
    }
    return sum;
}

/////////////////////////////////////////////////////////////////////////////////////


/*  ---------- Energy Calculation For Methane (CH4)  ---------- */
/*  could be potentially extended for any molecula with known:
*   1. number of degenerative vibrational modes (numberModes);
*   2. number of vibrational energy levels for each mode (numberVibrLvlByMode);
*   3. degrees of degeneraticy for each mode (dByMode).
*/

double OneTempApproxMultiModes::calcEnergy(macroParam& point)
{
    double UTrRot = getTrRotEnegry(point, 0);
    double UVibr = getVibrEnergy(point, 0);
    double E = point.density * (UTrRot + UVibr) + 0.5 * pow(point.velocity, 2) * point.density;
    return E;
}

double OneTempApproxMultiModes::getTrRotEnegry(macroParam& point, size_t component)
{
    double n = Nav * point.density / point.mixture.components[component].molarMass;
    double Utr = 3. / 2. * kB * point.temp * n / point.density;
    double Urot = kB * point.temp / point.mixture.components[component].mass;
    return Utr + Urot;
}

double OneTempApproxMultiModes::getVibrEnergy(macroParam& point, size_t component)
{
    double Uvibr = avgVibrEnergy(point, component) / point.mixture.mass(component);
    return Uvibr;
}

double OneTempApproxMultiModes::avgVibrEnergy(macroParam& point, size_t component)
{
    // fixed formulas from problem statement need to be added
    
    MixtureComponent molecula = point.mixture.components[component];
    double e_1000 = hc * molecula.omega_eByMode[0];// dimensions of omega were wrong!!
    double e_0100 = hc * molecula.omega_eByMode[1];
    double e_0010 = hc * molecula.omega_eByMode[2];
    double e_0001 = hc * molecula.omega_eByMode[3];
    double e_0000 = hc * (molecula.omega_eByMode[0] * molecula.dByMode[0] / 2. + molecula.omega_eByMode[1] * molecula.dByMode[1] / 2.
                          + molecula.omega_eByMode[2] * molecula.dByMode[2] / 2. + molecula.omega_eByMode[3] * molecula.dByMode[3] / 2.);

    double Z = Zvibr(point, component);

    double sum = 0;

    for (const auto& inds : point.mixture.components[component].possibleVibrInds) {
     
        double s = (inds[1] + 1) * (inds[2] + 1) * (inds[2] + 2) * (inds[3] + 1) * (inds[3] + 2) / 4.;

        double e_0 = inds[0] * e_1000 + inds[1] * e_0100 + inds[2] * e_0010 + inds[3] * e_0001;

        sum += s * (e_0 + e_0000) / (Z) * exp(-e_0 / (kB * point.temp));
    }
    return sum;
}

/* // don't use yet
double OneTempApproxMultiModes::vibrEnergyLvl(int lvl1, int lvl2, int lvl3, int lvl4, macroParam& point, size_t component)
{
    MixtureComponent molecula = point.mixture.components[component];

    double result = hc * (
        molecula.omega_eByMode[0] * (lvl1 + molecula.dByMode[0] / 2.) +
        molecula.omega_eByMode[1] * (lvl2 + molecula.dByMode[1] / 2.) +
        molecula.omega_eByMode[2] * (lvl3 + molecula.dByMode[2] / 2.) +
        molecula.omega_eByMode[3] * (lvl4 + molecula.dByMode[3] / 2.)
        );

    return result;
}*/


double OneTempApproxMultiModes::Zvibr(macroParam& point, size_t component)
{
    MixtureComponent molecula = point.mixture.components[component];
    double e_1000 = hc * molecula.omega_eByMode[0];
    double e_0100 = hc * molecula.omega_eByMode[1];
    double e_0010 = hc * molecula.omega_eByMode[2];
    double e_0001 = hc * molecula.omega_eByMode[3];

    double sum = 0;
    for (const auto& inds : point.mixture.components[component].possibleVibrInds) {
        double s = (inds[1] + 1) * (inds[2] + 1) * (inds[2] + 2) * (inds[3] + 1) * (inds[3] + 2) / 4.;
        double e_0 = inds[0] * e_1000 + inds[1] * e_0100 + inds[2] * e_0010 + inds[3] * e_0001;
        sum += s * exp(-e_0 / (kB * point.temp));
    }
    return sum;
}

double OneTempApproxMultiModes::getGamma(macroParam& point)
{
    /* function to calculate adiabatic index gamma = Cp/Cv */
    // mass, U=energy - v^2/2
    size_t component = 0; // we consider one-component methane gas
    double Cv_tr = 3.0 / 2 * kB / point.mixture.mass(component);
    double Cv_rot = kB / point.mixture.mass(component);

    double h = getEntalp(point, component);
    double nu_c = clight * point.mixture.components[component].omega_e;
    double Cv_vibr = (kB / point.mixture.mass(component)) * pow(h * nu_c / (kB * point.temp), 2) * exp(-h * nu_c / (kB * point.temp)); // (kB / point.mixture.mass(component)) * pow(h * nu_c / (kB * point.temp), 2) * exp(h * nu_c / (kB * point.temp)) / pow(exp(h * nu_c / (kB * point.temp)) - 1, 2);

    double Cv = Cv_tr + Cv_rot + Cv_vibr;

    double gamma = (UniversalGasConstant / point.mixture.molarMass(component) + Cv) / Cv;
    return gamma;
}

double OneTempApproxMultiModes::getEntalp(macroParam& point, size_t component)
{
    /* function to calculate specific entalpy h = p/rho + U */
    double UTrRot = getTrRotEnegry(point, 0);
    double UVibr = getVibrEnergy(point, 0);
    double res = 5. * kB * point.temp / (2. * point.mixture.mass(component)) + (UTrRot + UVibr);
    return res;
}
