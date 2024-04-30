#ifndef _MIXTURE
#define _MIXTURE
#include <vector>
#include <string>
#include "global.h"

struct MixtureComponent
{
	double molarMass; // молярная масса
	double epsilonDevK; // параметр в потенциале (для рассчёта вязкости)
	double Zinf; // параметр для модели Парка
	double mass;
	double sigma; // диаметр в метрах
	double omega_e; // спектроскопическая постоянная
	int numberVibrLvl; // число колебательных уровней
	int numberAtoms;

	double D_diss; // энергия диссоциации молекулы
	int numberOfModes; // число вырожденных мод
	std::vector<double> omega_eByMode; // спектроскопические постоянные для каждой вырожденной моды
	std::vector<int> numberVibrLvlByMode; // количество колебательных уровней для каждой вырожденной моды
	std::vector<int> dByMode; // степень вырожденности каждой моды
	std::vector<std::vector<int>> possibleVibrInds; // набор возможных состояний молекулы по колебательным уровням мод

	//... какие-то другие параметры компонент
	std::string name;                    //название компоненты
};

struct Mixture
{
	Mixture() { NumberOfComponents = 0; };
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
	double Zinf(size_t i);
};
#endif
