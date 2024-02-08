#pragma once
#include "global.h"
#include "macroparam.h"

struct BorderCondition
{
    virtual void updatePoints(vector<macroParam> &points) = 0;

    virtual double get_dyc_dy(){return 0;}; // затычка для более серьёзных условий
};

struct BorderConditionCouette : public BorderCondition
{
    void updatePoints(vector<macroParam> &points);
    void setWallParameters(double up_velocity_ , double down_velocity_ , double up_temp_ , double down_temp_)
        {up_velocity = up_velocity_ ; down_velocity = down_velocity_; up_temp = up_temp_; down_temp = down_temp_;};
    double get_dyc_dy(){return 0;};
protected:
    double up_velocity , down_velocity, up_temp , down_temp;
};

struct BorderConditionPersonal : public BorderConditionCouette
{
    void updatePoints(vector<macroParam> &points);
    double get_dyc_dy(){return 0;};
};

struct BorderConditionSoda : public BorderCondition
{
    void updatePoints(vector<macroParam> &points){return;};
    double get_dyc_dy(){return 0;};
};


struct BorderConditionShockwave : public BorderCondition
{
    void updatePoints(vector<macroParam>& points);
    void setWallParameters(double up_velocity_, double up_temp_, double up_pressure_, double up_density_)
    {
        up_velocity = up_velocity_; up_temp = up_temp_; up_pressure = up_pressure_; up_density = up_density_;
    };
    double get_dyc_dy() { return 0; };
protected:
    double up_velocity, up_temp, up_pressure, up_density;
};