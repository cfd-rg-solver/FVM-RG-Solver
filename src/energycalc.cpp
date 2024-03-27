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
    double E = point.density * (UTrRot + UVibr + 0.5 * pow(point.velocity, 2));
    /*
    if (point.temp == 300.) {
        std::cout << "-------------------------" << std::endl;
        std::cout << "T=" << point.temp << std::endl;
        std::cout << "p=" << point.pressure << std::endl;
        std::cout << "rho=" << point.density << std::endl;
        std::cout << "v=" << point.velocity << std::endl;
        double Ctr = 3. / 2. * kB / point.mixture.mass(0);
        double Crot = 3. / 2. * kB / point.mixture.mass(0);
        double Cvibr = getCvibr(point, 0);
        double Cv = Ctr + Crot + Cvibr;
        std::cout << "Cv=" << Cv << std::endl;
        std::cout << "zeta=" << getBulkViscosity(point,0) << std::endl;
        std::cout << "-------------------------" << std::endl;
    }*/
    return E;
}

double OneTempApproxMultiModes::getTrRotEnegry(macroParam& point, size_t component)
{
    double n = Nav * point.density / point.mixture.components[component].molarMass;
    double Utr = 3. / 2. * kB * point.temp * n / point.density;
    double Urot = 3. / 2. * kB * point.temp / point.mixture.components[component].mass; // because methane has 3 degrees of freedom
    return Utr + Urot;
}

double OneTempApproxMultiModes::getVibrEnergy(macroParam& point, size_t component)
{
    double Uvibr = avgVibrEnergy(point, component) / point.mixture.mass(component);
    return Uvibr;
}

double OneTempApproxMultiModes::avgVibrEnergy(macroParam& point, size_t component)
{
    MixtureComponent molecula = point.mixture.components[component];
    double e_1000 = hc * molecula.omega_eByMode[0];
    double e_0100 = hc * molecula.omega_eByMode[1];
    double e_0010 = hc * molecula.omega_eByMode[2];
    double e_0001 = hc * molecula.omega_eByMode[3];
    double e_0000 = hc * (
        molecula.omega_eByMode[0] * molecula.dByMode[0] / 2.
        + molecula.omega_eByMode[1] * molecula.dByMode[1] / 2.
        + molecula.omega_eByMode[2] * molecula.dByMode[2] / 2.
        + molecula.omega_eByMode[3] * molecula.dByMode[3] / 2.
        );

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

double OneTempApproxMultiModes::getEntalp(macroParam& point, size_t component)
{
    /* function to calculate specific entalpy h = p/rho + U */

    double UTrRot = getTrRotEnegry(point, 0);
    double UVibr = getVibrEnergy(point, 0);
    double res = 5. * kB * point.temp / (2. * point.mixture.mass(component)) + (UTrRot + UVibr);
    return res;
}


double OneTempApproxMultiModes::ZvibrDiff(macroParam& point, size_t component)
{
    /* function to calculate analytical partial derivative of Zvibr by T */

    MixtureComponent molecula = point.mixture.components[component];
    double e_1000 = hc * molecula.omega_eByMode[0];
    double e_0100 = hc * molecula.omega_eByMode[1];
    double e_0010 = hc * molecula.omega_eByMode[2];
    double e_0001 = hc * molecula.omega_eByMode[3];

    double sum = 0;
    for (const auto& inds : point.mixture.components[component].possibleVibrInds) {
        double s = (inds[1] + 1) * (inds[2] + 1) * (inds[2] + 2) * (inds[3] + 1) * (inds[3] + 2) / 4.;
        double e_0 = inds[0] * e_1000 + inds[1] * e_0100 + inds[2] * e_0010 + inds[3] * e_0001;
        sum += s * (e_0 / (kB * pow(point.temp, 2))) * exp(-e_0 / (kB * point.temp));
    }

    return sum;
}

double OneTempApproxMultiModes::getVibrEnergyDiff(macroParam& point, size_t component)
{
    /* function to calculate analytical partial derivative of Uvibr by T */

    MixtureComponent molecula = point.mixture.components[component];
    double e_1000 = hc * molecula.omega_eByMode[0];
    double e_0100 = hc * molecula.omega_eByMode[1];
    double e_0010 = hc * molecula.omega_eByMode[2];
    double e_0001 = hc * molecula.omega_eByMode[3];
    double e_0000 = hc * (
        molecula.omega_eByMode[0] * molecula.dByMode[0] / 2. 
        + molecula.omega_eByMode[1] * molecula.dByMode[1] / 2.
        + molecula.omega_eByMode[2] * molecula.dByMode[2] / 2. 
        + molecula.omega_eByMode[3] * molecula.dByMode[3] / 2.
        );

    double Z = Zvibr(point, component);
    double Zdiff = ZvibrDiff(point, component);

    double sum = 0;

    for (const auto& inds : point.mixture.components[component].possibleVibrInds) {

        double s = (inds[1] + 1) * (inds[2] + 1) * (inds[2] + 2) * (inds[3] + 1) * (inds[3] + 2) / 4.;

        double e_0 = inds[0] * e_1000 + inds[1] * e_0100 + inds[2] * e_0010 + inds[3] * e_0001;

        double derivativeZ = Zdiff / (-pow(Z, 2));
        double derivetiveExp = (e_0 / (kB * pow(point.temp,2))) / Z;
        sum += s * (e_0 + e_0000) * exp(-e_0 / (kB * point.temp)) * (derivativeZ + derivetiveExp);
    }
    return sum;
}

double OneTempApproxMultiModes::getCvibr(macroParam& point, size_t component)
{
    /* function to calculate Cv_vibr = (dU_vibr/dT)_V */

    double Cvibr = getVibrEnergyDiff(point, component) / point.mixture.mass(component); // partial derivative of Uvibr by T
    return Cvibr;
}

double OneTempApproxMultiModes::getBulkViscosity(macroParam& point, size_t component)
{
    /* function to calculate bulk viscosity zeta = kB*T/beta_int * (c_int/c_V)^2 */

    double Ctr = 3. / 2. * kB / point.mixture.mass(component);
    double Crot = 3. / 2. * kB / point.mixture.mass(component); // because methane has 3 degrees of freedom
    double Cvibr = getCvibr(point, component);
    double Cv = Ctr + Crot + Cvibr;
    double Cint = Cv - Ctr;

    double tau_rot = 1.1 * 10e-9; // atm.s ??? TODO
    double tau_vibr = 1.9 * 10e-6; // atm.s ??? TODO ok for 300K => 300^(-1/3)~0.149 ???
    double beta_int = point.pressure * UniversalGasConstant * pow(Cint/Cv,2) / (Crot/tau_rot + Cvibr/tau_vibr);

    double zeta = kB * point.temp * pow(Cint / Cv, 2) / (beta_int);
    return zeta;
}

double OneTempApproxMultiModes::getGamma(macroParam& point)
{
    /* function to calculate adiabatic index gamma = Cp/Cv */

    size_t component = 0; // we consider one-component methane gas, its fixed in order to use getGamma in systemOfEquation
    double Cv_tr = 3.0 / 2 * kB / point.mixture.mass(component);
    double Cv_rot = kB / point.mixture.mass(component);
    double Cv_vibr = getCvibr(point, component);
    double Cv = Cv_tr + Cv_rot + Cv_vibr;

    double gamma = (UniversalGasConstant / point.mixture.molarMass(component) + Cv) / Cv;
    return gamma;
}
