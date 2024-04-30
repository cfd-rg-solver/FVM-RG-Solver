#pragma once

//#include "abstractsolver.h"


struct HLLCSolver : public AbstractSolver
{
	HLLCSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_) :AbstractSolver(mixture_, startParam_, solParam_) {};

	// запускает процесс решения задачи
	void solve();


protected:

	//подготавливает размеры всех векторов
	void prepareVectors();

	// Расчет вектора потоков во всех ячейках
	void computeF();

	// Расчет релаксационных членов
	void computeR();

	// Расчет потоков на стыках ячеек методом HLLE
	void computeHlleF();

	// Расчет потоков на стыках ячеек методом HLLC
	void computeHllcF();

	// Расчет потоков на стыках ячеек методом HLL
	void computeHllF();

	void updateU();

	// обновлеяет вектор макропараметров с помощью U
	void updatePoints();


	// Значения потока на границах ячеек по методу HLLC
	Matrix  fluxF2, fluxF2_normal, fluxF3;
	vector<Matrix> fluxF1;

	//записывать ли данные в файл ?
	bool isWriteData = false;
};
