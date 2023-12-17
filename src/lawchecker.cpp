#include "lawchecker.h"
#include <iostream>
void LawChecker::rememberU(vector<double> &U, int t, double time)
{
    if(t == t1)
    {
        time1 = time;
        U1 = U;
        return;
    }
    if(t == t2)
    {
        time2 = time;
        U2 = U;
        checkLaw();
        return;
    }
    else
        return;
}

void LawChecker::rememberF(vector<double> &F_, vector<double> &Fv_, int t, double time)
{
    if(t == t1)
    {
        Fv = Fv_;
        F = F_;
    }
    return;
}

bool LawChecker::checkLaw()
{
    if(U1.size()*U2.size()*F.size() == 0)
        return false;

    double S1 = 0;
    double S2 = 0;
    for(int i = y1; i <= y2; i++)
    {
        S1 += U1[i];
        S2 += U2[i];
    }
    double deltaT = time2 - time1;
    double d = (S2 - S1) * dh + deltaT * (F[y2]+Fv[y2] - F[y1] - Fv[y1]);
    std::cout<<std::endl<< "DIFFERENCE "<< d << std::endl <<std::endl;
    return true;

}
