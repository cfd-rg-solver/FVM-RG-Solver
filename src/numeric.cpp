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
