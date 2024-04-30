#pragma once
#include "global.h"
#include "macroparam.h"

class Observer
{
public:
	Observer(double precision_) :precision(precision_) {};

	// задаёт через какое количество итераций наблюдатель должен проверять соотвествие
	void setPeriodicity(int period);

	// запоминает параметры на текущем шаге
	void remember(vector<macroParam> param);

	// возвращает true, если нужно продолжать рассчёт, false если параметры уже изменяются незначительно
	bool checkDifference(vector<macroParam> param);

	// возвращает periodicity
	int getPeriodicity();
private:
	int periodicity = 1;
	vector<macroParam> previousParam;
	double precision;
};
