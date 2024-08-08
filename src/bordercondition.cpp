#include "bordercondition.h"

void BorderConditionCouette::updatePoints(vector<macroParam>& points)
{
	bool presEq = 1;
	// bool slipBC = 1;
	size_t N = points.size();
	Mixture mixture = points[1].mixture;

	points[0].mixture = mixture;
	if (presEq)
		points[0].pressure = points[1].pressure;
	else
		points[0].density = points[1].density;
	points[0].densityArray = points[1].densityArray;
	points[0].fractionArray = points[1].fractionArray;

	points[0].velocity_tau = -points[1].velocity_tau + 2. * down_velocity;
	points[0].velocity_normal = -points[1].velocity_normal;

	points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau), 2) + pow(fabs(points[0].velocity_normal), 2));
	points[0].temp = -points[1].temp + 2. * down_temp;

	// дополнительные рассчитываемые величины

	if (!presEq)
		points[0].pressure = points[0].density / mixture.molarMass(points[0].fractionArray) * UniversalGasConstant * points[0].temp;
	else
		points[0].density = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (UniversalGasConstant * points[0].temp);

	//solParam.NumCell-1
	points[N - 1].mixture = mixture;
	if (presEq)
		points[N - 1].pressure = points[N - 2].pressure;
	else
		points[N - 1].density = points[N - 2].density;
	points[N - 1].densityArray = points[N - 2].densityArray;
	points[N - 1].fractionArray = points[N - 2].fractionArray;
	points[N - 1].velocity_tau = -points[N - 2].velocity_tau + 2. * up_velocity;
	points[N - 1].velocity_normal = -points[N - 2].velocity_normal;
	points[N - 1].velocity = sqrt(pow(points[N - 1].velocity_tau, 2) + pow(points[N - 1].velocity_normal, 2));
	points[N - 1].temp = -points[N - 2].temp + 2. * up_temp;
	// дополнительные рассчитываемые величины
	auto y_c = points[N - 1].fractionArray;
	if (!presEq)
		points[N - 1].pressure = points[N - 1].density / mixture.molarMass(y_c) * UniversalGasConstant * points[N - 1].temp;
	else
		points[N - 1].density = points[N - 1].pressure * mixture.molarMass(y_c) / (UniversalGasConstant * points[N - 1].temp);
}

void BorderConditionPersonal::updatePoints(vector<macroParam>& points)
{

	bool presEq = 1;
	size_t N = points.size();
	Mixture mixture = points[1].mixture;
	//0
	points[0].mixture = mixture;

	points[0].pressure = points[1].pressure;

	points[0].density = points[1].density;
	points[0].densityArray = points[1].densityArray;
	points[0].fractionArray = points[1].fractionArray;

	points[0].velocity_tau = points[1].velocity_tau;
	//    points[0].velocity_normal = points[1].velocity_normal;
	points[0].velocity_normal = points[1].velocity_normal;


	points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau), 2) + pow(fabs(points[0].velocity_normal), 2));

	points[0].temp = points[0].pressure / (points[0].density * UniversalGasConstant) * mixture.molarMass(points[0].fractionArray);

	// дополнительные рассчитываемые величины

	//solParam.NumCell-1
	points[N - 1].mixture = mixture;
	if (presEq)
		points[N - 1].pressure = points[N - 2].pressure;
	else
		points[N - 1].density = points[N - 2].density;
	points[N - 1].densityArray = points[N - 2].densityArray;
	points[N - 1].fractionArray = points[N - 2].fractionArray;
	points[N - 1].velocity_tau = -points[N - 2].velocity_tau + 2. * up_velocity;
	points[N - 1].velocity_normal = -points[N - 2].velocity_normal;
	points[N - 1].velocity = sqrt(pow(points[N - 1].velocity_tau, 2) + pow(points[N - 1].velocity_normal, 2));
	points[N - 1].temp = -points[N - 2].temp + 2. * up_temp;
	// дополнительные рассчитываемые величины
	auto y_c = points[N - 1].fractionArray;
	if (!presEq)
		points[N - 1].pressure = points[N - 1].density / mixture.molarMass(y_c) * UniversalGasConstant * points[N - 1].temp;
	else
		points[N - 1].density = points[N - 1].pressure * mixture.molarMass(y_c) / (UniversalGasConstant * points[N - 1].temp);
}

void BorderConditionShockwave::updatePoints(vector<macroParam>& points)
{
	int BCtype = 1;
	size_t N = points.size();

	points[0].velocity_normal = points[1].velocity_normal;
	points[0].velocity_tau = 0;
	points[0].velocity = points[0].velocity_normal;

	Mixture mixture = points[1].mixture;

	if (BCtype == 0) {
		points[0].velocity_normal = 0;
		points[0].velocity_tau = left_velocity;
		points[0].velocity = left_velocity;

		points[0].density = left_density;
		points[0].densityArray[0] = left_density; // !!! only for 1 component gas
		points[0].fractionArray = points[1].fractionArray;

		points[0].temp = left_temp;
		points[0].pressure = points[0].density * UniversalGasConstant * points[0].temp / mixture.molarMass(points[0].fractionArray);

		points[N - 1].velocity_normal = 0;
		points[N - 1].velocity_tau = right_velocity;
		points[N - 1].velocity = right_velocity;

		points[N - 1].density = right_density;
		points[N - 1].densityArray[0] = right_density; // !!! only for 1 component gas
		points[N - 1].fractionArray = points[N - 2].fractionArray;

		points[N - 1].temp = right_temp;
		points[N - 1].pressure = points[N - 1].density * UniversalGasConstant * points[N - 1].temp / mixture.molarMass(points[N - 1].fractionArray);
	}
	else if (BCtype == 1)
	{
		points[0].velocity_normal = 0;
		points[0].velocity_tau = points[1].velocity_tau;
		points[0].velocity = points[1].velocity_tau;

		points[0].densityArray = points[1].densityArray;
		points[0].fractionArray = points[1].fractionArray;

		points[0].pressure = points[1].pressure;
		points[0].temp = points[1].temp; // из уравнения состояния ид газа
		points[0].density = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (points[0].temp * UniversalGasConstant);

		points[N - 1].velocity_normal = 0;
		points[N - 1].velocity_tau = points[N - 1].velocity_tau;
		points[N - 1].velocity = points[N - 1].velocity_tau;

		points[N - 1].densityArray = points[N - 2].densityArray;
		points[N - 1].fractionArray = points[N - 2].fractionArray;

		points[N - 1].pressure = points[N - 2].pressure;
		points[N - 1].temp = points[N - 2].temp; // из уравнения состояния ид газа
		points[N - 1].density = points[N - 1].pressure * mixture.molarMass(points[N - 1].fractionArray) / (points[N - 1].temp * UniversalGasConstant);
	}
	else if (BCtype == 2) 
	{
		points[0].velocity_normal = 0;
		points[0].velocity_tau = left_velocity;
		points[0].velocity = left_velocity;

		points[0].densityArray = points[1].densityArray;
		points[0].fractionArray = points[1].fractionArray;

		points[0].pressure = points[1].pressure;
		points[0].density = points[1].density;
		points[0].temp = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (points[0].density * UniversalGasConstant); // из уравнения состояния ид газа

		points[N - 1].velocity_normal = 0;
		points[N - 1].velocity_tau = right_velocity;
		points[N - 1].velocity = right_velocity;

		points[N - 1].densityArray = points[N - 2].densityArray;
		points[N - 1].fractionArray = points[N - 2].fractionArray;

		points[N - 1].pressure = points[N - 2].pressure;
		points[N - 1].density = points[N - 2].density;
		points[N - 1].temp = points[N - 1].pressure * mixture.molarMass(points[N - 1].fractionArray) / (points[N - 1].density * UniversalGasConstant); // из уравнения состояния ид газа
	}
	points[0].gamma = energyCalculator->getGamma(points[0]);
	points[N - 1].gamma = energyCalculator->getGamma(points[N - 1]);
}

void BorderConditionSoda::updatePoints(vector<macroParam>& points)
{
	size_t N = points.size();
	Mixture mixture = points[1].mixture;
	points[0].velocity_normal = 0;
	points[0].velocity_tau = points[1].velocity_tau;
	points[0].velocity = points[0].velocity_tau;

	points[0].densityArray = points[1].densityArray;
	points[0].fractionArray = points[1].fractionArray;

	points[0].pressure = points[1].pressure;
	points[0].density = points[1].density;
	points[0].temp = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (points[0].density * UniversalGasConstant); // из уравнения состояния ид газа

	points[N - 1].velocity_normal = 0;
	points[0].velocity_tau = points[N - 2].velocity_tau;
	points[N - 1].velocity = points[N - 1].velocity_tau;

	points[N - 1].densityArray = points[N - 2].densityArray;
	points[N - 1].fractionArray = points[N - 2].fractionArray;

	points[N - 1].pressure = points[N - 2].pressure;
	points[N - 1].density = points[N - 2].density;
}
