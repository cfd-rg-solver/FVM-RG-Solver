#pragma once
#include <cmath>
#include "macroparam.h"
double newtonMethod( double (*func) (macroParam), double (*dfunc) (macroParam), macroParam p0, double rp )
{
  double t1  = p0.temp - (func(p0) - rp) / dfunc(p0); // первое приближение
  double eps = 0.0001;
  while (fabs(t1 - p0.temp) > eps)
  {
    p0.temp = t1;
    t1 = p0.temp - (func(p0) - rp) / dfunc(p0); // последующие приближения
  }
  return t1;
}
