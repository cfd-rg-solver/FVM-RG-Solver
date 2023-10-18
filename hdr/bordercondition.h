#pragma once
#include "global.h"
struct BorderCondition
{
    double up_velocity , down_velocity = 0., up_temp , down_temp;

    double get_dyc_dy(); //затычка для более серьёзных условий
};

struct BorderConditionSoda  // временное решение
{
    double leftDensity, leftPressure, leftVelocity, rightDensity, rightPressure, rightVelocity;
};
