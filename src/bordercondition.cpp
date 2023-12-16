#include "bordercondition.h"

void BorderConditionCouette::updatePoints(vector<macroParam> &points)
{

    bool presEq = 1;
    bool slipBC = 1; // TO DO
    Mixture mixture = points[0].mixture;
    int N = points.size();
    //0
    if(presEq)
        points[0].pressure = points[1].pressure;
    else
        points[0].density = points[1].density;
    points[0].densityArray = points[1].densityArray;
    points[0].fractionArray = points[1].fractionArray;
    points[0].velocity_tau = - points[1].velocity_tau + 2.* down_velocity;
    points[0].velocity_normal = - points[1].velocity_normal;
    points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
    points[0].temp = - points[1].temp +  2.* down_temp;
    // дополнительные рассчитываемые величины
    if(!presEq)
        points[0].pressure = points[0].density / mixture.molarMass(points[0].fractionArray) * UniversalGasConstant * points[0].temp;
    else
        points[0].density = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (UniversalGasConstant * points[0].temp);
    points[0].soundSpeed = sqrt(gamma*points[0].pressure/points[0].density);


    //solParam.NumCell-1
    if(presEq)
        points[N-1].pressure = points[N-2].pressure;
    else
        points[N-1].density = points[N-2].density;
    points[N-1].densityArray = points[N-2].densityArray;
    points[N-1].fractionArray = points[N-2].fractionArray;
    points[N-1].velocity_tau = - points[N-2].velocity_tau + 2.* up_velocity;
    points[N-1].velocity_normal = - points[N-2].velocity_normal;
    points[N-1].velocity = sqrt(pow(points[N-2].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
    points[N-1].temp = - points[N-2].temp +  2.* up_temp;
    // дополнительные рассчитываемые величины
    auto y_c = points[N - 1].fractionArray;
    if(!presEq)
        points[N - 1].pressure = points[N - 1].density / mixture.molarMass(y_c) * UniversalGasConstant * points[N-1].temp;
    else
        points[N - 1].density = points[N - 1].pressure * mixture.molarMass(y_c) / (UniversalGasConstant * points[N-1].temp);
    points[N - 1].soundSpeed = sqrt(gamma*points[N-1].pressure/points[N-1].density);

}

