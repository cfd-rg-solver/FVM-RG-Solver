#pragma once
#include "global.h"
#include "macroparam.h"

struct BorderCondition
{
    virtual void updatePoints(vector<macroParam> points) = 0;

    virtual double get_dyc_dy(){return 0;}; //затычка для более серьёзных условий

    void setGamma(double gamma_){gamma = gamma_;};

    double up_velocity , down_velocity = 0., up_temp , down_temp;
    double gamma;
};

struct BorderConditionCouette : public BorderCondition
{
    void updatePoints(vector<macroParam> points);

    double get_dyc_dy(); //затычка для более серьёзных условий
};

struct BorderConditionSoda : public BorderCondition
{
    void updatePoints(vector<macroParam> points){return;};
};
