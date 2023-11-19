#ifndef _MIXTURE
#define _MIXTURE
#include <vector>
#include <string>
#include "global.h"
struct MixtureComponent
{
    double density;                 //плотность компоненты
    double molarMass; // молярная масса
    double epsilonDevK; // параметр в потенциале (для рассчёта вязкости) (хз пока чё это такое)
    double mass;
    double sigma; // диаметр в метрах
    double omega_e; // спектроскопическая постоянная
    int numberVibrLvl; // число колебательных уровней
    int numberAtoms;
    //... какие-то другие параметры компонент
    std::string name;                    //название компоненты

    double avgVibrEnergyDiff(double temp);
    double avgVibrEnergy(double temp);
    double vibrEnergyLvl(int lvl);
    double ZvibrDiff(double temp);
    double Zvibr(double temp);


};
struct Mixture
{
    Mixture(){NumberOfComponents = 0;};
    Mixture(std::vector<MixtureComponent> components_);

    std::vector<MixtureComponent> components;
    int NumberOfComponents;
    //... какие-то другие параметры смеси
    //пример

    double molarMass();
    double molarMass(size_t i);
    double mass(size_t i);
    double sigma(size_t i);
    double epsilonDevK(size_t i);
};
#endif
