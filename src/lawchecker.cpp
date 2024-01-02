#include "lawchecker.h"
#include <iostream>
#include <algorithm>
void LawChecker::rememberU(vector<double> &U,vector<double> &rho, int t, double time)
{
    if(t == t1)
    {
        time1 = time;
        U1 = U;
        rho1 = rho;
        return;
    }
    if(t == t2)
    {
        time2 = time;
        U2 = U;
        rho2 = rho;
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

    double N = U1.size();
    double I1 = 0;
    double I2 = 0;
    for(int k = y1+1; k < y2 - 1; k+=2)
    {
        I1+= U1[k-1] + 4*U1[k] + U1[k+1];
        I2+= U2[k-1] + 4*U2[k] + U2[k+1];
    }
    I1 *= dh / 3.;
    I2 *= dh / 3.;

    vector<double> dU1(N-1  ),dU2(N-1  );
    vector<double> ddU1(N-2 ),ddU2(N-2 );
    vector<double> dddU1(N-3),dddU2(N-3);
    for(int i = 0; i < dU1.size(); i++)
    {
        dU1[i] = (U1[i+1] - U1[i]) / dh;
        dU2[i] = (U2[i+1] - U2[i]) / dh;
    }
    for(int i = 0; i < ddU1.size(); i++)
    {
        ddU1[i] = (dU1[i+1] - dU1[i]) / dh;
        ddU2[i] = (dU2[i+1] - dU2[i]) / dh;
    }
    for(int i = 0; i < dddU1.size(); i++)
    {
        dddU1[i] = (ddU1[i+1] - ddU1[i]) / dh;
        dddU2[i] = (ddU2[i+1] - ddU2[i]) / dh;
    }
    auto mU1 = max_element(std::begin(dddU1), std::end(dddU1));
    auto mU2 = max_element(std::begin(dddU2), std::end(dddU2));
    double coeff = dh*(N - 1)/288. * pow(dh,3);
    double Ef1 = *mU1 * coeff;
    double Ef2 = *mU2 * coeff;

    double deltaT = time2 - time1;
    double d = (I2 - I1) + deltaT * (F[y2]+Fv[y2] - F[y1] - Fv[y1]);

    std::cout<<std::endl<< "DIFFERENCE "<< d << std::endl;
    std::cout<<" U|t2 - U|t1 = " <<(I2 - I1)<<std::endl;
    std::cout<<" dF = "<< deltaT * (F[y2]+Fv[y2] - F[y1] - Fv[y1])<<std::endl;
    std::cout<< "Ef1 = "<< Ef1 <<std::endl;
    std::cout<< "Ef2 = "<< Ef2 <<std::endl;
    std::cout<< "dT = "<<deltaT<<std::endl;
    vector<double> dRho(N);
    for(int i = 0; i < dRho.size();i++)
    {
        dRho[i] = abs(rho2[i] - rho1[i]);
    }
    std::cout<< "max(deltaRho) = "<< *max_element(std::begin(dRho), std::end(dRho))<<std::endl<<std::endl;
    return true;

}
