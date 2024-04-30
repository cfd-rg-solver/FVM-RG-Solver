#include <functional>

#include "abstractsolver.h"

AbstractSolver::AbstractSolver(Mixture mixture_, solverParams solParam_, SystemOfEquationType type, RiemannSolverType riemannType)
{
	isContinue = 0;
	mixture = mixture_;
	solParam = solParam_;
	delta_h = 0;

	if (mixture.NumberOfComponents == 1)
	{
		auto* tmp = new CoeffSolver1Comp1Temp();
		coeffSolver = tmp;
	}
	else
	{
		auto* tmp = new CoeffSolver2Comp1Temp();
		coeffSolver = tmp;
	}

	system = getSystemOfEquation(type);
	riemannSolver = getRiemannSolver(riemannType);
	riemannSolver->solParam = solParam_;

	system->setCoeffSolver(coeffSolver);
	system->setMixture(mixture);
	system->setNumberOfCells(solParam.NumCell);
	system->setSolverParams(solParam);
	system->setEqSolver(eqSolver);
}

void AbstractSolver::setWriter(DataWriter* writer_)
{
	writer = writer_;
	isWriteData = true;
}

void AbstractSolver::setObserver(Observer* obs)
{
	observer = obs;
	watcherIteration = observer->getPeriodicity();
	isObserverWatching = true;
}

void AbstractSolver::setStartDistribution(StartCondition* startDist)
{
	prepareVectorSizes();
	startDist->setStartDistribution(points);
	mixture = points[0].mixture;
	system->prepareSolving(points);
}


void AbstractSolver::writePoints(double i)
{
	if (isWriteData)
		writer->writeData(points, i);
}

void AbstractSolver::setDelta_h(double dh)
{
	delta_h = dh;
}

void AbstractSolver::setBorderConditions(BorderCondition* border_)
{
	border = border_;
	system->setBorderCondition(border);
}

void AbstractSolver::correctData()
{
	for (size_t i = 0; i < points.size(); i++)
	{
		points[i].mixture = mixture;
		points[i].soundSpeed = sqrt(solParam.Gamma * points[i].pressure / points[i].density);
	}
	return;
}

SystemOfEquation* AbstractSolver::getSystemOfEquation(SystemOfEquationType type)
{
	switch (type)
	{
	case SystemOfEquationType::couette2:
	{
		auto* tmp = new Couette2();
		return tmp;
	}
	case SystemOfEquationType::couette2Alt:
	{
		auto* tmp = new Couette2Alt();
		return tmp;
	}
	case SystemOfEquationType::couette2AltBinary:
	{
		auto* tmp = new Couette2AltBinary();
		return tmp;
	}
	case SystemOfEquationType::soda:
	{
		auto* tmp = new Soda();
		return tmp;
	}
	case SystemOfEquationType::shockwave1:
	{
		auto* tmp = new Shockwave1();
		return tmp;
	}
	case SystemOfEquationType::shockwave2:
	{
		auto* tmp = new Shockwave2();
		return tmp;
	}
	}
	return nullptr;
}

RiemannSolver* AbstractSolver::getRiemannSolver(RiemannSolverType type)
{
	switch (type)
	{
	case RiemannSolverType::HLLCSolver:
	{
		return new struct HLLCSolver();
	}
	case RiemannSolverType::HLLESolver:
	{
		return new struct HLLESolver();
	}
	case RiemannSolverType::HLLSimple:
	{
		return new struct HLLSimple();
	}

	case RiemannSolverType::ExacRiemanSolver:
	{
		return new struct ExacRiemanSolver();
	}
	}
	return nullptr;
}

void AbstractSolver::setEnergyCalculator(EnergyCalc* energyCalculator_)
{
	energyCalculator = energyCalculator_;
	system->setEnergyCalculator(energyCalculator);
};

void AbstractSolver::prepareVectorSizes()
{
	points.resize(solParam.NumCell);
	for (size_t i = 0; i < points.size(); i++)
	{
		points[i].densityArray.resize(mixture.NumberOfComponents);
		points[i].fractionArray.resize(mixture.NumberOfComponents);
	}
	system->prepareVectorSizes();
}


void AbstractSolver::setDt()
{
	double max = riemannSolver->maxSignalVelocity;  // это нужно чтобы правильно подобрать временной шаг, чтобы соблюдался критерий КФЛ
	double dt;
	if (max != 0)
		dt = solParam.CFL * delta_h / max;
	else
		dt = 1e-7;// 0.00001;
	timeSolvind.push_back(dt);
	return;
}

void AbstractSolver::updatePoints()
{
	auto size = points.size() - 1;
	// #pragma omp parallel for schedule(static)
	for (int i = 0; i < size + 1; i++)
		//    for(int i = 1; i < size; i++)
	{
		points[i].density = system->getDensity(i);
		points[i].velocity_tau = system->getVelocityTau(i);
		points[i].velocity_normal = system->getVelocityNormal(i);
		points[i].velocity = system->getVelocity(i);
		points[i].pressure = system->getPressure(i);

		for (size_t j = 0; j < mixture.NumberOfComponents; j++)
		{
			points[i].densityArray[j] = system->getDensity(i, j); // тут косяк TODO
			points[i].fractionArray[j] = points[i].densityArray[j] / points[i].density;
		}
		points[i].soundSpeed = sqrt(solParam.Gamma * points[i].pressure / points[i].density);
		points[i].temp = system->getTemp(i);
		points[i].gamma = energyCalculator->getGamma(points[i]);
	}
	border->updatePoints(points);
}


bool AbstractSolver::observerCheck(size_t currentIteration)
{
	if (currentIteration % watcherIteration == 0)
	{
		observer->remember(points);
		return true;
	}
	if (currentIteration % watcherIteration == 1)
	{
		return observer->checkDifference(points);
	}
	return true;
}

bool AbstractSolver::lawCheck(size_t currentIteration)
{
	size_t checkerIteration = 10000;
	if (currentIteration % checkerIteration == 0)
	{
		lawChecker.collectData(system->U, system->F, system->Fv, timeSolvind.last(), 0);
		return true;
	}
	if (currentIteration % checkerIteration == 1)
	{
		lawChecker.collectData(system->U, system->F, system->Fv, timeSolvind.last(), 1);
		cout << "Energy error = " << lawChecker.checkEnergyLaw() << endl;
		return true;
	}
	return true;
}