#include "mixture.h"
#include "global.h"

Mixture::Mixture(std::vector<MixtureComponent> components_)
{
    components = components_;
    NumberOfComponents = components.size();
}

double Mixture::molarMass()
{
    double sum = 0;
    for(size_t i = 0; i < NumberOfComponents; i++)
        sum += components[i].molarMass;
    return sum;
}

double Mixture::molarMass(std::vector<double> y_c)
{
    double sum = 0;
    for(size_t j = 0; j < y_c.size(); j++)
    {
        double M_c = molarMass(j);
        sum += y_c[j] / M_c;
    }
    double M = 1/sum;
    return M;
}


double Mixture::molarMass(size_t i)
{
    return components[i].molarMass;
}

double Mixture::mass(size_t i)
{
    return components[i].mass;
}

double Mixture::sigma(size_t i)
{
    return components[i].sigma;
}

double Mixture::epsilonDevK(size_t i)
{
    return components[i].epsilonDevK;
}


double MixtureComponent::avgVibrEnergyDiff(double temp)
{
    if(numberAtoms == 1)
        return 0;

    double Z = Zvibr(temp);
    double Zdiff = ZvibrDiff(temp);

    double tmp = 0;
    double tmpDiff = 0;
    double sum = 0;
    for(size_t i = 0; i < numberVibrLvl; i++)
    {
        double eps_ic = vibrEnergyLvl(i);
        tmp = eps_ic * exp (-eps_ic / ( kB * temp));
        tmpDiff = pow(eps_ic,2) * exp (-eps_ic / ( kB * temp)) / ( kB * pow(temp,2));
        sum += tmpDiff * Z - tmp * Zdiff;
    }

    double res = sum / (pow( Z , 2));
    return res;
}

double MixtureComponent::avgVibrEnergy(double temp)
{
    if(numberAtoms == 1)
        return 0;

    double sum = 0;
    double Z = Zvibr(temp);
    for(size_t i = 0; i < numberVibrLvl; i++)
    {
        sum += vibrEnergyLvl(i) / Z * exp (-vibrEnergyLvl(i) / ( kB * temp));
    }
    return sum;
}

double MixtureComponent::vibrEnergyLvl(int lvl)
{
    return hc * omega_e * (lvl + 1./2.);
}

double MixtureComponent::ZvibrDiff(double temp)
{
    if(numberAtoms == 1)
        return 0;

    double sum = 0;
    double s_i = 1; // только для O2
    for(size_t i = 0; i < numberVibrLvl; i++)
    {
        double eps_ic = vibrEnergyLvl(i);
        sum += s_i * eps_ic * exp (-eps_ic / ( kB * temp)) / ( kB * pow(temp , 2.));
    }
    return sum;
}

double MixtureComponent::Zvibr(double temp)
{
    if(numberAtoms == 1)
        return 0;

    double sum = 0;
    double s_i = 1; // только для O2
    for(size_t i = 0; i < numberVibrLvl; i++)
    {
        sum += exp (-vibrEnergyLvl(i) / ( kB * temp));
    }
    return sum;
}
