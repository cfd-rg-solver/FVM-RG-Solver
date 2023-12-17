#pragma once
#include <vector>

using std::vector;
struct LawChecker
{
    LawChecker(int y1_,int y2_,int t1_,int t2_, double dh_):y1(y1_), y2(y2_), t1(t1_), t2(t2_), dh(dh_){};

    int y1 = 0, y2 = 0;
    int t1 = 0, t2 = 0;

    double dh;

    double time1,time2;

    void rememberU(vector<double> &U, int t, double time);
    void rememberF(vector<double> &F_, vector<double> &Fv_, int t, double time);
    bool checkLaw();
private:
    vector<double> U1,U2,F,Fv;

};
