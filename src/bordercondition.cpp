#include "bordercondition.h"

void BorderConditionCouette::updatePoints(vector<macroParam> &points)
{
    bool presEq = 1;
    size_t N = points.size();
    Mixture mixture = points[1].mixture;
    if(presEq)
    {
        //0
        points[0].mixture = mixture;
        points[0].pressure = points[1].pressure;
        points[0].fractionArray = points[1].fractionArray;
        points[0].velocity_tau = -points[1].velocity_tau + 2.* down_velocity;
        points[0].velocity_normal = -points[1].velocity_normal;
        points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
        points[0].temp = -points[1].temp +  2. * down_temp;
        points[0].density = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (UniversalGasConstant * points[0].temp);

        for(int i = 0; i < points[0].mixture.components.size(); i++)
            points[0].densityArray[i] = points[0].density * points[0].fractionArray[i];


        //solParam.NumCell-1
        points[N-1].mixture = mixture;
        points[N-1].pressure = points[N-2].pressure;
        points[N-1].fractionArray = points[N-2].fractionArray;
        points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* up_velocity;
        points[N-1].velocity_normal = -points[N-2].velocity_normal;
        points[N-1].velocity = sqrt(pow(points[N-1].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
        points[N-1].temp = -points[N-2].temp +  2.* up_temp;
        points[N-1].density = points[N-1].pressure * mixture.molarMass(points[N-1].fractionArray) / (UniversalGasConstant * points[N-1].temp);

        for(int i = 0; i < points[N-1].mixture.components.size(); i++)
            points[N-1].densityArray[i] = points[N-1].density * points[N-1].fractionArray[i];
    }
    else
    {
        //TODO
    }
    return;
}

void BorderConditionPersonal::updatePoints(vector<macroParam> &points)
{

    bool presEq = 1;
    size_t N = points.size();
    Mixture mixture = points[1].mixture;
    //0
    points[0].mixture = mixture;

    points[0].pressure = points[1].pressure;

    points[0].density = points[1].density;
    points[0].densityArray = points[1].densityArray;
    points[0].fractionArray = points[1].fractionArray;

    points[0].velocity_tau = points[1].velocity_tau;
    //    points[0].velocity_normal = points[1].velocity_normal;
    points[0].velocity_normal = points[1].velocity_normal;


    points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));

    points[0].temp = points[0].pressure / ( points[0].density * UniversalGasConstant) * mixture.molarMass(points[0].fractionArray);

    // дополнительные рассчитываемые величины

    //solParam.NumCell-1
    points[N-1].mixture = mixture;
    if(presEq)
        points[N-1].pressure = points[N-2].pressure;
    else
        points[N-1].density = points[N-2].density;
    points[N-1].densityArray = points[N-2].densityArray;
    points[N-1].fractionArray = points[N-2].fractionArray;
    points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* up_velocity;
    points[N-1].velocity_normal = -points[N-2].velocity_normal;
    points[N-1].velocity = sqrt(pow(points[N-1].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
    points[N-1].temp = -points[N-2].temp +  2.* up_temp;
    // дополнительные рассчитываемые величины
    auto y_c = points[N-1].fractionArray;
    if(!presEq)
        points[N-1].pressure = points[N-1].density / mixture.molarMass(y_c) * UniversalGasConstant * points[N-1].temp;
    else
        points[N-1].density = points[N-1].pressure * mixture.molarMass(y_c) / (UniversalGasConstant * points[N-1].temp);
}

void BorderConditionShockwave::updatePoints(vector<macroParam>& points)
{
    bool BCtype = 1;
    size_t N = points.size();
    Mixture mixture = points[1].mixture;

    if (BCtype) {
        points[0].velocity_normal = 0;
        points[0].velocity_tau = left_velocity;
        points[0].velocity = left_velocity;

        points[0].densityArray = points[1].densityArray;
        points[0].fractionArray = points[1].fractionArray;

        points[0].pressure = points[1].pressure;
        points[0].density = points[1].density;
        points[0].temp = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (points[0].density * UniversalGasConstant); // из уравнения состояния ид газа

        points[N-1].velocity_normal = 0;
        points[N-1].velocity_tau = right_velocity;
        points[N-1].velocity = right_velocity;

        points[N-1].densityArray = points[N-2].densityArray;
        points[N-1].fractionArray = points[N-2].fractionArray;

        points[N-1].pressure = points[N-2].pressure;
        points[N-1].density = points[N-2].density;
        points[N-1].temp = points[N-1].pressure * mixture.molarMass(points[N-1].fractionArray) / (points[N-1].density * UniversalGasConstant); // из уравнения состояния ид газа
    }

    // ! another option (check if correct)
    else { // not changed
        points[0].velocity_normal = left_velocity;
        points[0].velocity_tau = 0;
        points[0].velocity = points[0].velocity_normal;

        points[0].densityArray = points[1].densityArray;
        points[0].fractionArray = points[1].fractionArray;

        points[0].density = left_density;
        points[0].temp = left_temp;
        points[0].pressure = points[0].density * UniversalGasConstant * points[0].temp / mixture.molarMass(points[0].fractionArray);

        points[N-1].velocity_normal = right_velocity;
        points[N-1].velocity_tau = 0;
        points[N-1].velocity = points[N-1].velocity_normal;

        points[N-1].densityArray = points[N-2].densityArray;
        points[N-1].fractionArray = points[N-2].fractionArray;

        points[N-1].density = right_density;
        points[N-1].temp = right_temp;
        points[N-1].pressure = points[N-1].density * UniversalGasConstant * points[N-1].temp / mixture.molarMass(points[N-1].fractionArray) ;
    }

}

void BorderConditionSoda::updatePoints(vector<macroParam> &points)
{
    size_t N = points.size();
    Mixture mixture = points[1].mixture;
    points[0].velocity_normal = 0;
    points[0].velocity_tau = points[1].velocity_tau;
    points[0].velocity = points[0].velocity_tau;

    points[0].densityArray = points[1].densityArray;
    points[0].fractionArray = points[1].fractionArray;

    points[0].pressure = points[1].pressure;
    points[0].density = points[1].density;

    points[N-1].velocity_normal = 0;
    points[0].velocity_tau = points[N-2].velocity_tau;
    points[N-1].velocity = points[N-1].velocity_tau;

    points[N-1].densityArray = points[N-2].densityArray;
    points[N-1].fractionArray = points[N-2].fractionArray;

    points[N-1].pressure = points[N-2].pressure;
    points[N-1].density = points[N-2].density;
}

void BorderConditionCouetteSlip::updatePointsStart(vector<macroParam> &points)
{
    bool presEq = 1;
    size_t N = points.size();
    Mixture mixture = points[1].mixture;
    if(presEq)
    {
        //0
        points[0].mixture = mixture;
        points[0].pressure = points[1].pressure;
        points[0].fractionArray = points[1].fractionArray;
        points[0].velocity_tau = -points[1].velocity_tau + 2.* down_velocity;
        points[0].velocity_normal = -points[1].velocity_normal;
        points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
        points[0].temp = -points[1].temp +  2. * down_temp;
        points[0].density = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (UniversalGasConstant * points[0].temp);

        for(int i = 0; i < points[0].mixture.components.size(); i++)
            points[0].densityArray[i] = points[0].density * points[0].fractionArray[i];


        //solParam.NumCell-1
        points[N-1].mixture = mixture;
        points[N-1].pressure = points[N-2].pressure;
        points[N-1].fractionArray = points[N-2].fractionArray;
        points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* up_velocity;
        points[N-1].velocity_normal = -points[N-2].velocity_normal;
        points[N-1].velocity = sqrt(pow(points[N-1].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
        points[N-1].temp = -points[N-2].temp +  2.* up_temp;
        points[N-1].density = points[N-1].pressure * mixture.molarMass(points[N-1].fractionArray) / (UniversalGasConstant * points[N-1].temp);

        for(int i = 0; i < points[N-1].mixture.components.size(); i++)
            points[N-1].densityArray[i] = points[N-1].density * points[N-1].fractionArray[i];
    }
}

void BorderConditionCouetteSlip::updatePoints(vector<macroParam> &points)
{
    // NOW CORRECT ONLY FOR SINGLE COMPONENT GAS !!!!!!!!!!!!
    bool presEq = 1;
    size_t N = points.size();
    Mixture mixture = points[1].mixture;
    if(presEq)
    {
        //0
        points[0].mixture = mixture;
        points[0].pressure = points[1].pressure;
        points[0].fractionArray = points[1].fractionArray;
        points[0].density = points[0].pressure * mixture.molarMass(points[0].fractionArray) / (UniversalGasConstant * points[0].temp);
        for(int i = 0; i < points[0].mixture.components.size(); i++)
            points[0].densityArray[i] = points[0].density * points[0].fractionArray[i];

        double velocityHalf;
        double temperatureHalf;

        velocityHalf = calcVelocityHalf(points[0], points[1], 0, down_velocity);
        temperatureHalf = calcTempHalf(points[0], points[1], 0, down_velocity, down_temp, velocityHalf);

        points[0].velocity_tau = -points[1].velocity_tau + 2.* velocityHalf;
        points[0].velocity_normal = -points[1].velocity_normal;
        points[0].velocity = sqrt(pow(fabs(points[0].velocity_tau),2) + pow(fabs(points[0].velocity_normal),2));
        points[0].temp = -points[1].temp + temperatureHalf;



        //solParam.NumCell-1
        points[N-1].mixture = mixture;
        points[N-1].pressure = points[N-2].pressure;
        points[N-1].fractionArray = points[N-2].fractionArray;
        points[N-1].density = points[N-1].pressure * mixture.molarMass(points[N-1].fractionArray) / (UniversalGasConstant * points[N-1].temp);
        for(int i = 0; i < points[N-1].mixture.components.size(); i++)
            points[N-1].densityArray[i] = points[N-1].density * points[N-1].fractionArray[i];

        velocityHalf = calcVelocityHalf(points[N-1], points[N-2], 0, up_velocity);
        temperatureHalf = calcTempHalf(points[N-1], points[N-2], 0, up_velocity, up_temp, velocityHalf);

        points[N-1].velocity_tau = -points[N-2].velocity_tau + 2.* velocityHalf;
        points[N-1].velocity_normal = -points[N-2].velocity_normal;
        points[N-1].velocity = sqrt(pow(points[N-1].velocity_tau,2) + pow(points[N-1].velocity_normal,2));
        points[N-1].temp = -points[N-2].temp +  2.* temperatureHalf;
    }
    else
    {
        //TODO
    }
}

double BorderConditionCouetteSlip::calcVelocityHalf(macroParam p0, macroParam p1, size_t component, double wallVelocity)
{
    // p0 - ghost cell, p1 - real cell
    int i = component;
    double sigma = p0.mixture.sigma(i);
    double rhoHalf = p0.densityArray[i];
    double m = p0.mixture.mass(i);
    double M = p0.mixture.molarMass(i);
    macroParam point(p0.mixture);
    point.temp = interp1(p0.temp, p1.temp);
    for(size_t j = 0; j < point.mixture.NumberOfComponents; j++)
    {
        point.fractionArray[j] = interp1(p0.densityArray[j], p1.densityArray[j]);
    }
    double mu = coeffSolver->shareViscositySimple(point);
    double T_last = point.temp;
    double mult = sqrt(2 * M_PI/ m * kB * T_last) * (2 - sigma) * kB *  M / (sigma * UniversalGasConstant * rhoHalf);
    double numerator = -mult * p1.velocity_tau / delta_h + wallVelocity;
    double denominator = 1 - mult / delta_h;
    return numerator / denominator;
}

double BorderConditionCouetteSlip::calcTempHalf(macroParam p0, macroParam p1, size_t component, double wallVelocity, double wallTemperature, double velocityHalf)
{
    // p0 - ghost cell, p1 - real cell
    int i = component;
    double sigma = p0.mixture.sigma(i);
    double rhoHalf = p0.densityArray[i];
    double m = p0.mixture.mass(i);
    double M = p0.mixture.molarMass(i);
    macroParam point(p0.mixture);
    point.temp = interp1(p0.temp, p1.temp);
    for(size_t j = 0; j < point.mixture.NumberOfComponents; j++)
    {
        point.fractionArray[j] = interp1(p0.densityArray[j], p1.densityArray[j]);
    }
    double T_last = point.temp;
    double mult = (2 - sigma) / (2 * sigma) * sqrt((M_PI * m) / (2 * kB * T_last)) * (M * coeffSolver->lambda(point)) / (UniversalGasConstant * rhoHalf);
    double add = m * pow((velocityHalf - wallVelocity),2) / (4 * kB);
    double numerator = -mult * p1.temp / delta_h + wallTemperature + add ;
    double denominator = 1 - mult / delta_h;
    return numerator / denominator;
}

double BorderConditionCouetteSlip::interp1(double value1, double value2)
{
    return (value1 + value2) / 2;
}
