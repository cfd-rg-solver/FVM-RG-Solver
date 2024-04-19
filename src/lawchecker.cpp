#include "lawchecker.h"
#include "numeric.h"

void LawChecker::collectData(std::vector<std::vector<double> > U_, std::vector<std::vector<double> > F_, std::vector<std::vector<double> > Fv_, double time_, bool isNew)
{
    if (isNew)
    {
        next.makeScreen(U_, F_, Fv_, time_);
    }
    else
    {
        prev.makeScreen(U_, F_, Fv_, time_);
    }
}

double LawChecker::checkEnergyLaw()
{
    IntegralCalculatorTrap calculator;
    size_t energyEqIdx = prev.U.size() - 1;
    double E2 = calculator.integrate(next.U[energyEqIdx], dh);
    double E1 = calculator.integrate(prev.U[energyEqIdx], dh);

    double integralPart = E2 - E1;

    double dt = next.time - prev.time;
    double resF1 = prev.resF[energyEqIdx].front();// at y1
    double resF2 = prev.resF[energyEqIdx].back(); // at y2

    double derrivativePart = -dt * (resF2 - resF1);

    double error = integralPart - derrivativePart;
    return error;
}

void FlowScreen::makeScreen(std::vector<std::vector<double> > U_, std::vector<std::vector<double> > F_, std::vector<std::vector<double> > Fv_, double time_)
{
    U = U_;
    F = F_;
    Fv = Fv_;
    time = time_;

    resF = F;

    truncateParameter(U);
    truncateParameter(F);
    truncateParameter(resF);

    if (Fv.empty())
        return;
    truncateParameter(Fv);

    for (int i = 0; i < U.size(); i++)
    {
        for (int j = 0; j < U[0].size(); j++)
            resF[i][j] += Fv[i][j];
    }

    return;
}

void FlowScreen::truncateParameter(std::vector<std::vector<double> >& vec)
{
    for (int i = 0; i < vec.size(); i++)
    {
        vec[i].pop_back();
        vec.erase(vec.begin());
    }
    return;
}