#include <functional>

#include "abstractsolver.h"

AbstractSolver::AbstractSolver(Mixture mixture_, solverParams solParam_, SystemOfEquationType type, RiemannSolverType riemannType)
{
    isContinue = 0;
    mixture = mixture_;
    solParam =solParam_;
    delta_h = 0;

    if(mixture.NumberOfComponents == 1)
    {
        auto *tmp = new CoeffSolver1Comp1Temp();
        coeffSolver = tmp;
    }
    else
    {
        auto *tmp = new CoeffSolver2Comp1Temp();
        coeffSolver = tmp;
    }

    system = getSystemOfEquation(type);
    riemannSolver = getRiemannSolver(riemannType);
    riemannSolver->solParam = solParam_;

    system->setBorderCondition(border);
    system->setCoeffSolver(coeffSolver);
    system->setMixture(mixture);
    system->setNumberOfCells(solParam.NumCell);
    system->setSolverParams(solParam);
    system->setEqSolver(eqSolver);
}


void AbstractSolver::setBorderConditions(double up_velocity_, double up_temp_, double down_temp_)
{
    border = new BorderConditionCouette();
    border->setGamma(solParam.Gamma);
    border->up_velocity = up_velocity_;
    border->down_velocity = 0.;
    border->up_temp =  up_temp_;
    border->down_temp = down_temp_;
    return;
}

void AbstractSolver::setBorderConditions()
{
    border = new BorderConditionSoda();
    border->setGamma(solParam.Gamma);
    return;
}

void AbstractSolver::setWriter(DataWriter *writer_)
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

void AbstractSolver::setStartDistribution(vector<macroParam> start)
{
    prepareVectorSizes();
    mixture = start[1].mixture;
    points = start;
    system->prepareSolving(points);
}

void AbstractSolver::setStartDistribution(macroParam start)
{
    prepareVectorSizes();
    mixture = start.mixture;
    for(size_t i = 1; i < points.size()-1; i++)
    {
        points[i].mixture = mixture;
        points[i].temp = start.temp;
        points[i].fractionArray =  start.fractionArray;
        points[i].density = start.density;

        points[i].pressure = points[i].density * UniversalGasConstant * start.temp/mixture.molarMass(start.fractionArray);
        //points[i].density = startParam.pressure * mixture.molarMass()/(UniversalGasConstant * startParam.temp);
        points[i].densityArray =  start.densityArray;
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
        points[i].velocity_tau = start.velocity_tau;
        points[i].velocity_normal = 0;
        points[i].velocity = fabs(points[i].velocity_tau);
    }
    // для points[0] и points[solParam.NumCell-1] (!важно что идёт после цикла!)
    useBorder();
    system->prepareSolving(points);
}

void AbstractSolver::setStartDistribution(macroParam left, macroParam right)
{
    mixture = left.mixture;
    prepareVectorSizes();
    for(size_t i = 0; i < points.size()/2; i++)
    {
        points[i].mixture = mixture;
        points[i].pressure = left.pressure;
        points[i].density  = left.density;
        points[i].velocity_tau = left.velocity;
        points[i].velocity_normal = 0;
        points[i].velocity = left.velocity;
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
    }
    for(size_t i = points.size()/2; i < points.size(); i++)
    {
        points[i].mixture = mixture;
        points[i].pressure = right.pressure;
        points[i].density  = right.density;
        points[i].velocity_tau = right.velocity;
        points[i].velocity_normal = 0;
        points[i].velocity = right.velocity;
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
    }
    system->prepareSolving(points);
}

void AbstractSolver::writePoints(double i)
{
    if(isWriteData)
        writer->writeData(points,i);
}

void AbstractSolver::setDelta_h(double dh)
{
    delta_h = dh;
}

void AbstractSolver::correctData()
{
    for(size_t i = 0; i < points.size(); i++)
    {
        points[i].mixture = mixture;
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
    }
    return;
}

SystemOfEquation *AbstractSolver::getSystemOfEquation(SystemOfEquationType type)
{
    switch(type)
    {
        case SystemOfEquationType::couette2:
        {
            auto *tmp = new Couette2();
            return tmp;
        }
        case SystemOfEquationType::couette2Alt:
        {
            auto *tmp = new Couette2Alt();
            return tmp;
        }
        case SystemOfEquationType::couette2AltBinary:
        {
            auto *tmp = new Couette2AltBinary();
            return tmp;
        }
        case SystemOfEquationType::soda:
        {
            auto *tmp = new Soda();
            return tmp;
        }
    }
    return nullptr;
}

RiemannSolver *AbstractSolver::getRiemannSolver(RiemannSolverType type)
{
    switch(type)
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

void AbstractSolver::prepareVectorSizes()
{
    points.resize(solParam.NumCell);
    for(size_t i = 0; i < points.size(); i++)
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
    if(max!=0)
        dt = solParam.CFL*delta_h/max;
    else
        dt = 0.00001;
    timeSolvind.push_back(dt);
    return;
}

void AbstractSolver::updatePoints()
{
    auto size = points.size()-1;
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size + 1; i++)
//    for(int i = 1; i < size; i++)
    {
        points[i].velocity_tau = system->getVelocityTau(i);
        points[i].velocity_normal = system->getVelocityNormal(i);
        points[i].velocity = system->getVelocity(i);
        points[i].density = system->getDensity(i);
        points[i].pressure = system->getPressure(i);

        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        {
            points[i].densityArray[j] =  system->getDensity(i,j); // тут косяк TODO
            points[i].fractionArray[j] = points[i].densityArray[j] / points[i].density;
        }
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
        points[i].temp = system->getTemp(i);
    }
    useBorder();
}

//void AbstractSolver::useBorder()
//{
//    //0
//    points[0].mixture = mixture;
//    points[0].pressure = points[1].pressure;
//    points[0].densityArray =points[1].densityArray;
//    points[0].fractionArray =points[1].fractionArray;
//    points[0].velocity_tau = -points[1].velocity_tau + 2.*border.down_velocity;
//    points[0].velocity_normal = -points[1].velocity_normal;
//    points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
//    points[0].temp = -points[1].temp +  2.*border.down_temp;
//    // дополнительные рассчитываемые величины
//    points[0].density = points[0].pressure * mixture.molarMass() /UniversalGasConstant / points[0].temp;
//    points[0].soundSpeed = sqrt(solParam.Gamma*points[0].pressure/points[0].density);


//    //solParam.NumCell-1
//    points[solParam.NumCell-1].mixture = mixture;
//    points[solParam.NumCell-1].pressure = points[solParam.NumCell-2].pressure;
//    points[solParam.NumCell-1].densityArray = points[solParam.NumCell-2].densityArray;
//    points[solParam.NumCell-1].fractionArray = points[solParam.NumCell-2].fractionArray;
//    points[solParam.NumCell-1].velocity_tau = -points[solParam.NumCell-2].velocity_tau + 2.*border.up_velocity;
//    points[solParam.NumCell-1].velocity_normal = -points[solParam.NumCell-2].velocity_normal;
//    points[solParam.NumCell-1].velocity = sqrt(pow(points[solParam.NumCell-1].velocity_tau,2) + pow(points[solParam.NumCell-1].velocity_normal,2));
//    points[solParam.NumCell-1].temp = -points[solParam.NumCell-2].temp +  2.*border.up_temp;
//    // дополнительные рассчитываемые величины
//    points[solParam.NumCell-1].density = points[solParam.NumCell-1].pressure * mixture.molarMass() /UniversalGasConstant / points[solParam.NumCell-1].temp;
//    points[solParam.NumCell-1].soundSpeed = sqrt(solParam.Gamma*points[solParam.NumCell-1].pressure/points[solParam.NumCell-1].density);
//}

void AbstractSolver::useBorder()
{
   border->updatePoints(points);
}


bool AbstractSolver::observerCheck(size_t currentIteration)
{
    if(currentIteration%watcherIteration == 0)
    {
        observer->remember(points);
        return true;
    }
    if(currentIteration%watcherIteration == 1)
    {
        return observer->checkDifference(points);
    }
    return true;
}

