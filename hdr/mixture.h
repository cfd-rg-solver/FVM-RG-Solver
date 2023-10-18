#ifndef _MIXTURE
#define _MIXTURE
#include <vector>
#include <string>

struct MixtureComponent
{
    double density;                 //плотность компоненты
    double molarMass; // молярная масса
    double epsilonDevK; // параметр в потенциале (для рассчёта вязкости) (хз пока чё это такое)
    double mass;
    double sigma; // диаметр в метрах
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

    double molarMass();
    double molarMass(size_t i);
    double mass(size_t i);
    double sigma(size_t i);
    double epsilonDevK(size_t i);
    double getEffDiff(size_t i); // пока временное решение
    double getEntalp(size_t i);
};
#endif
