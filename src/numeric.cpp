#include "numeric.h"

double Newton::solveEq(EnergyCalc *energy , macroParam point, double rightPart)
{
    macroParam p2 = point;
    macroParam p3 = point;
    double dh = 0.0001;
    p2.temp += dh;
    p3.temp -= dh;

    auto f = energy->getEnergyFunc();
    double f2 = f(p2);
    double f3 = f(p3);
    double df23 = (f2 - f3) / (2 * dh);

    double t1  = point.temp - (f(point) - rightPart) / df23; // первое приближение
    double eps = dh / 10;
    while (fabs(t1 - point.temp) > eps)
    {
        point.temp = t1;
        p2.temp = t1 + 0.00001;
        p3.temp = t1 - 0.00001;
        f2 = f(p2);
        f3 = f(p3);
        df23 = (f2 - f3) / 0.00002;

        t1 = point.temp - (f(point) - rightPart) / df23; // последующие приближения
    }
    return t1;
}

double IntegralCalculatorTrap::integrate(std::vector<double> x, std::vector<double> y)
{
    double sum = 0;
    for (size_t i = 1; i < y.size(); i++)
    {
        sum += (y[i - 1] + y[i]) / 2. * (x[i] - x[i - 1]);
    }
    return sum;
}
double IntegralCalculatorTrap::integrate(std::vector<double> y, double dh)
{
    double sum = 0;
    for (size_t i = 1; i < y.size(); i++)
    {
        sum += (y[i - 1] + y[i]) / 2. * dh;
    }
    return sum;
}

double DerrivativeCalculator::derrivativeCenter(std::vector<double> y, double dh, int idx)
{
    if (idx == 0 || idx == y.size() - 1)
        return 0;
    return (y[idx + 1] - y[idx - 1]) / (2 * dh);
}

double DerrivativeCalculator::derrivativeLeft(std::vector<double> y, double dh, int idx)
{
    if (idx == 0)
        return 0;
    return (y[idx] - y[idx - 1]) / (dh);
}

double DerrivativeCalculator::derrivativeRight(std::vector<double> y, double dh, int idx)
{
    if (idx == y.size() - 1)
        return 0;
    return (y[idx + 1] - y[idx]) / (dh);
}