#pragma once
#include <vector>

struct FlowScreen
{
    double time;
    std::vector<std::vector<double>> U, F, Fv, resF;
    void makeScreen(std::vector<std::vector<double>> U_, std::vector<std::vector<double>> F_, std::vector<std::vector<double>> Fv_, double time_);

private:
    void truncateParameter(std::vector<std::vector<double>>& vec);
};

struct LawChecker
{
    double dh;
    FlowScreen prev, next;
    void setDh(double dh_) { dh = dh_; };
    void collectData(std::vector<std::vector<double>> U_, std::vector<std::vector<double>> F_, std::vector<std::vector<double>> Fv_, double time_, bool isNew);
    double checkEnergyLaw();
};