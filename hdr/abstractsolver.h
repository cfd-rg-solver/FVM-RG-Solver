#pragma once

#include <vector>
#include <list>
#include <mutex>

#include "global.h"
#include "mixture.h"
#include "coeffsolver.h"
#include "bordercondition.h"
#include "startcondition.h"
#include "datawriter.h"
#include "observer.h"
#include "systemofequation.h"
#include "riemannsolver.h"
#include "numeric.h"
#include "lawchecker.h"

struct AbstractSolver
{
public:
	AbstractSolver(Mixture mixture_, solverParams solParam_, SystemOfEquationType type, RiemannSolverType riemannType);
	// запускает процесс решения задачи
	virtual void solve() = 0;
	//устанавливает размер ячейки
	void setDelta_h(double dh);

	// устанавливает некоторые граничные условия (TODO сделать более общую структуру)
	void setBorderConditions(BorderCondition* border_);

	// устанавливает записыватель и поднимает флаг записи
	void setWriter(DataWriter* writer_);

	// устанавливает наблюдателя, который будет следить за тем, чтобы рассчёт остановился, когда будет достигнута новая точность
	void setObserver(Observer* obs);

	// устанавливает начальное распределение
	void setStartDistribution(StartCondition* startDist);

	// устанавливает энергетические расчеты
	void setEnergyCalculator(EnergyCalc* energyCalculator_);

	SystemOfEquation* system;

	RiemannSolver* riemannSolver;

	Observer* observer;

	Mixture mixture;

	solverParams solParam;

	CoeffSolver* coeffSolver;

	LawChecker lawChecker;

	NonLinearEqSolver* eqSolver = new Newton();

	BorderCondition* border;

	DataWriter* writer;

	EnergyCalc* energyCalculator;

protected:

	//записывает текущие макропараметры points[] в папку с названием i
	void writePoints(double i);

	// функция для дозаполнения полей вектора points значениями скорости звука и параметрами смеси, которые не вытаскиваются из файлов расчёта
	virtual void correctData();

	SystemOfEquation* getSystemOfEquation(SystemOfEquationType type);
	RiemannSolver* getRiemannSolver(RiemannSolverType type);

	//подготавливает размеры всех векторов
	virtual void prepareVectorSizes();
	//    //подготавливает размеры всех векторов по заданным параметрам
	//    virtual void prepareVectorSizes();

	// устанавливает временной шаг
	void setDt();

	void updatePoints();

	bool lawCheck(size_t currentIteration);

	void UpdateBorderU();

	//функция для работы с наблюдателем, выдаёт true если нужно продолжать рассчёт
	bool observerCheck(size_t currentIteration);

	// хранит значение макропараметров в каждой ячейке
	vector<macroParam> points;

	//    //скорость звука в каждой ячейке
	//    Matrix sound_speed;

	// шаг сетки
	double delta_h = 0;

	// вектор в котором хранятся временные шаги
	Matrix timeSolvind;

	//записывать ли данные в файл ?
	bool isWriteData = false;

	//смотрит ли наблюдатель за установлением равновесия?
	bool isObserverWatching = false;

	// через какое количество итераций наблюдатель должен проверять соотвествие
	//(берутся две соседних итерации, а не те, что разделены watcherIteration итерациями)
	int watcherIteration = 1;

	// продолжается ли рассчёт? 0 - расчёт начинается с нуля, 1 - расчёт продолжается согласно данным из считанного файла
	bool isContinue = 0;
};
