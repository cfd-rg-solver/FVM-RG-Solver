#ifndef _MIXTURE
#define _MIXTURE
#include <vector>
#include <string>
#include "global.h"
struct MixtureComponent
{
    double molarMass; // молярная масса
    double epsilonDevK; // параметр в потенциале (для рассчёта вязкости)
    double mass;
    double sigma; // диаметр в метрах
    double omega_e; // спектроскопическая постоянная
    int numberVibrLvl; // число колебательных уровней
    int numberAtoms;
    //... какие-то другие параметры компонент
    std::string name;                    //название компоненты
};
struct Mixture
{
    Mixture(){NumberOfComponents = 0;};
    Mixture(std::vector<MixtureComponent> components_);

    std::vector<MixtureComponent> components;
    int NumberOfComponents;
    //... какие-то другие параметры смеси
    //пример

    double molarMass(); // неправильный рассчет в случае двух компонентного газа, использовать незя
    double molarMass(std::vector<double> y_c);
    double molarMass(size_t i);
    double mass(size_t i);
    double sigma(size_t i);
    double epsilonDevK(size_t i);
};
#endif
