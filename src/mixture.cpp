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

double Mixture::Zinf(size_t i)
{
    return components[i].Zinf;
}