#include"startcondition.h"

void UniformDistribution::setStartDistribution(vector<macroParam>& points)
{
	if (newSolving)
	{
		for (size_t i = 0; i < points.size(); i++)
		{
			points[i].mixture = example.mixture;
			points[i].temp = example.temp;
			points[i].fractionArray = example.fractionArray;
			points[i].density = example.density;
			points[i].pressure = points[i].density * UniversalGasConstant * example.temp / example.mixture.molarMass(example.fractionArray);
			//points[i].density = startParam.pressure * mixture.molarMass()/(UniversalGasConstant * startParam.temp);
			points[i].densityArray = example.densityArray;
			points[i].velocity_tau = example.velocity_tau;
			points[i].velocity_normal = example.velocity_normal;
			points[i].velocity = fabs(points[i].velocity_tau);
		}
	}
	else
	{
		points = exampleVec;
	}
	return;
}



void UniformDistributionBorder::setStartDistribution(vector<macroParam>& points)
{
	points[0].mixture = example.mixture;
	points[points.size() - 1].mixture = example.mixture;
	for (size_t i = 1; i < points.size() - 1; i++)
	{
		points[i].mixture = example.mixture;
		points[i].temp = example.temp;
		points[i].fractionArray = example.fractionArray;
		points[i].density = example.density;

		points[i].pressure = points[i].density * UniversalGasConstant * example.temp / example.mixture.molarMass(example.fractionArray);
		//points[i].density = startParam.pressure * mixture.molarMass()/(UniversalGasConstant * startParam.temp);
		points[i].densityArray = example.densityArray;
		points[i].velocity_tau = example.velocity_tau;
		points[i].velocity_normal = example.velocity_normal;
		points[i].velocity = fabs(points[i].velocity_tau);
	}
	// для points[0] и points[solParam.NumCell-1] (!важно что идёт после цикла!)
	borderCondition->updatePoints(points); // usingBorderCondition
	return;
}

void GapDistribution::setStartDistribution(vector<macroParam>& points)
{
	Mixture mixture = left.mixture;
	bool startType = 1;
	if (startType) {
		for (size_t i = 0; i < points.size() / 2 + 1; i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = left.pressure;
			points[i].density = left.density;
			points[i].fractionArray = left.fractionArray;
			points[i].densityArray = left.densityArray;
			points[i].velocity_tau = left.velocity_tau;
			points[i].velocity_normal = left.velocity_normal;
			points[i].velocity = left.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass() / points[i].density / UniversalGasConstant;
		}
		for (size_t i = points.size() / 2 + 1; i < points.size(); i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = right.pressure;
			points[i].density = right.density;
			points[i].fractionArray = right.fractionArray;
			points[i].densityArray = right.densityArray;
			points[i].velocity_tau = right.velocity_tau;
			points[i].velocity_normal = right.velocity_normal;
			points[i].velocity = right.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass() / points[i].density / UniversalGasConstant;
		}
	}
	else {
		for (size_t i = 0; i < 1; i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = left.pressure;
			points[i].density = left.density;
			points[i].fractionArray = left.fractionArray;
			points[i].densityArray = left.densityArray;
			points[i].velocity_tau = left.velocity_tau;
			points[i].velocity_normal = left.velocity_normal;
			points[i].velocity = left.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass() / points[i].density / UniversalGasConstant;
		}
		for (size_t i = 1; i < points.size(); i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = right.pressure;
			points[i].density = right.density;
			points[i].fractionArray = right.fractionArray;
			points[i].densityArray = right.densityArray;
			points[i].velocity_tau = right.velocity_tau;
			points[i].velocity_normal = right.velocity_normal;
			points[i].velocity = right.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass() / points[i].density / UniversalGasConstant;
		}
	}
}



void ShockwaveGapDistribution::setStartDistribution(vector<macroParam>& points)
{
	Mixture mixture = left.mixture;
	bool startType = 1;
	if (startType) {
		for (size_t i = 0; i < points.size() / 2 + 1; i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = left.pressure;
			points[i].density = left.density;
			points[i].fractionArray = left.fractionArray;
			points[i].densityArray = left.densityArray;
			points[i].velocity_tau = left.velocity_tau;
			points[i].velocity_normal = left.velocity_normal;
			points[i].velocity = left.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass() / points[i].density / UniversalGasConstant;
			points[i].gamma = energyCalculator->getGamma(points[i]);
		}
		for (size_t i = points.size() / 2 + 1; i < points.size(); i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = right.pressure;
			points[i].density = right.density;
			points[i].fractionArray = right.fractionArray;
			points[i].densityArray = right.densityArray;
			points[i].velocity_tau = right.velocity_tau;
			points[i].velocity_normal = right.velocity_normal;
			points[i].velocity = right.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass() / points[i].density / UniversalGasConstant;
			points[i].gamma = energyCalculator->getGamma(points[i]);
		}
	}
	else {
		for (size_t i = 0; i < 1; i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = left.pressure;
			points[i].density = left.density;
			points[i].fractionArray = left.fractionArray;
			points[i].densityArray = left.densityArray;
			points[i].velocity_tau = left.velocity_tau;
			points[i].velocity_normal = left.velocity_normal;
			points[i].velocity = left.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass() / points[i].density / UniversalGasConstant;
			points[i].gamma = energyCalculator->getGamma(points[i]);
		}
		for (size_t i = 1; i < points.size(); i++)
		{
			points[i].mixture = mixture;
			points[i].pressure = right.pressure;
			points[i].density = right.density;
			points[i].fractionArray = right.fractionArray;
			points[i].densityArray = right.densityArray;
			points[i].velocity_tau = right.velocity_tau;
			points[i].velocity_normal = right.velocity_normal;
			points[i].velocity = right.velocity;
			points[i].temp = points[i].pressure * points[i].mixture.molarMass() / points[i].density / UniversalGasConstant;
			points[i].gamma = energyCalculator->getGamma(points[i]);
		}
	}
}