#pragma once

#include <vector>
#include <cmath>
#include <string>

using namespace std;

const double Nav(6.02214129e23);
const double UniversalGasConstant = 8.3144598;
const double kB = 1.38064852e-23;    // J/K
const double kBE = 8.617e-5;       // const Bolz elVolt*K
static const double massaCO2 = 7.306e-26; // m_CO_2
static const double gasConst = kB / massaCO2;
static const double e100 = 2.757135054E-20;
static const double e010 = 1.32493446E-20;
static const double e001 = 4.66607366E-20;
static const double De = 8.83859e-19;

static const double hPlank = 6.62559e-34;
static const double clight = 2.99792458e8;
static const double hc = hPlank * clight; // hc = ww!!!
static const double o2[3] = { 1345.04e2, 667.25e2, 2361.71e2 }; // 1/m!Herzberg
static const double ox[7] = { -3.63e2, 3.44e2, -19.28e2, -0.635e2, -12.51e2, -12.56e2, 0.775e2 }; // 1/m!Herzberg ox11, ox12, ox13, ox22, ox23, ox33, oxll
static const double sigmaCO2CO2_mix2 = 3.763e-10; // collision diameter in m
static const double alpha1_CO2 = 17.5 / sigmaCO2CO2_mix2;
static const double a_SSH[3] = { 1. / 2., 8. / 11., 1. / 2. };
static const double a1_SSH[3] = { 8. / 11., 1. / 2., 3. / 11. };

static const double masRed_CO2_CO2 = massaCO2 / 2.;
static const double m[3] = { 1.46e-26, 1.338e-26, 1.46e-26 };
static const double ww = 1.60219e-21 / 8065.47; // if you need energy in J (not in meV) 1.60219e-19 / 8065.47
static const double ni[3] = { o2[0] * ww / hPlank, o2[1] * ww / hPlank, o2[2] * ww / hPlank };


const  double ZettaInf = 20.39;
extern double sigma;
extern double epsilonDevK;
extern double molMass;
extern double mass;

struct solverParams
{
	int NumCell = 0;    // Число расчтеных ячеек
	double t_fin = 0;    // Время выхода из решения
	double Ma = 0;    // Число маха
	double Gamma = 0;    // Показатель адиабаты
	double CFL = 0;    // Число Куранта
	double lambda = 0;    // Длина свободного пробега
	int lambdaSol = 0;    // Кол-во длин пробега для расчета
	int PlotIter = 0;    // Кол-во итераци, через которое отрисывыввается график
	int MaxIter = 10000000; // максимальное кол-во шагов по времени
	int typePlot = 0;    // Тип отображения на графике :
	// 0 - абс. давлени
	// 1 - абс. плотноть
	// 2 - абс. скорость
	// 3 - абс. температура
	int typeRightBorder = 1;//  Тип граничного условия справа
	int typeLeftBorder = 1; //  Тип граничного условия слева
};
class Matrix
{
private:
	std::vector<double> data;
public:
	Matrix(vector<double> val)
	{
		data = val;
	}
	Matrix()
	{

	}
	void clear()
	{
		data.clear();
	}
	Matrix(int len, double val = 0)
	{
		data = vector<double>(len, val);
	}
	vector<double>::iterator begin()
	{
		return data.begin();
	}
	vector<double>::iterator end()
	{
		return data.end();
	}
	int size()
	{
		return data.size();
	}
	operator vector<double>() const
	{
		return data;
	}
	double& operator [](int i)
	{
		static double elseValue;
		if (i < data.size())
			return data[i];
		else
			return elseValue;
	}
	double first()
	{
		return data.front(); // changed
	}
	double last()
	{
		return data.back(); // changed
	}
	void push_back(double val)
	{
		data.push_back(val);
	}
	void push_front(double val)
	{
		data.insert(data.begin(), val);  // changed
	}
	void removeLast()
	{
		data.pop_back();
	}
	void removeFirst()
	{
		data.erase(data.begin()); // changed
	}
	void resize(int i)
	{
		data.resize(i);
	}
	Matrix operator /(vector<double> div)
	{
		vector<double> output;
		if (data.size() != div.size())
			return data;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] / div[i]);
		return output;
	}
	Matrix operator *(vector<double> div)
	{
		vector<double> output;
		if (data.size() != div.size())
			return data;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] * div[i]);
		return output;
	}
	Matrix operator +(vector<double> div)
	{
		vector<double> output;
		if (data.size() != div.size())
			return data;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] + div[i]);
		return output;
	}
	Matrix operator -(vector<double> div)
	{
		vector<double> output;
		if (data.size() != div.size())
			return data;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] - div[i]);
		return output;
	}
	Matrix operator /(double div)
	{
		vector<double> output;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] / div);
		return output;
	}
	Matrix operator *(double div)
	{
		vector<double> output;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] * div);
		return output;
	}
	Matrix operator +(double div)
	{
		vector<double> output;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] + div);
		return output;
	}
	Matrix operator -(double div)
	{
		vector<double> output;
		for (int i = 0; i < data.size(); i++)
			output.push_back(data[i] - div);
		return output;
	}
	static Matrix POW(vector<double> div, double param)
	{
		vector<double> output;
		for (int i = 0; i < div.size(); i++)
			output.push_back(pow(div[i], param));
		return output;
	}
	static Matrix SQRT(vector<double> div)
	{
		vector<double> output;
		for (int i = 0; i < div.size(); i++)
			output.push_back(sqrt(div[i]));
		return output;
	}
	static Matrix REVERSE(vector<double> div)
	{
		vector<double> output;
		for (int i = 0; i < div.size(); i++)
			output.push_back(1.0 / div[i]);
		return output;
	}
};
