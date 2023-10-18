#include "godunovsolversoda.h"
#include <iostream>
#include <omp.h>
//#include <typeinfo.h>
void GodunovSolverSoda::solve()
{
    prepareSolving();
    writePoints(-1);
    double T = 0;
    for(size_t i  = 0; i < solParam.MaxIter; i++)
    {
        // Устанавливаем текущий временной шаг
        setDt();
        T += timeSolvind.last();

        system->computeF(points, delta_h);

        riemannSolver->computeFlux(system);
        //riemannSolver->computeFlux(system,timeSolvind.last(), delta_h);

        // Вычисляем вектор релаксационных членов
        //computeR();

        // Обновляем вектор U
        system->updateU(delta_h,timeSolvind.last());
        // Обновляем вектор макропараметров
        updatePoints();


        //записать данные, если это требуется
        //writePoints(T*1000000); // микросек
        if(i%1 == 0)
        {
            std::cout<<i<<" iteration"<<std::endl;
            writePoints(i); // микросек
        }

        //проверка точности
        if(isObserverWatching)
        {
            // то есть если проверка наблюдателя не пройдена, нужно прекратить рассчёт
            if(!observerCheck(i))
                break;
        }
    }
    writePoints(T*1000000);
}

void GodunovSolverSoda::prepareSolving()
{
    prepareVectorSizes();
    for(size_t i = 0; i < points.size()/2; i++)
    {
        points[i].mixture = mixture;
        points[i].pressure = borderSoda.leftPressure;
        points[i].density  = borderSoda.leftDensity;
        points[i].velocity_tau = borderSoda.leftVelocity;
        points[i].velocity_normal = 0;
        points[i].velocity = borderSoda.leftVelocity;
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
    }
    for(size_t i = points.size()/2; i < points.size(); i++)
    {
        points[i].mixture = mixture;
        points[i].pressure = borderSoda.rightPressure;
        points[i].density  = borderSoda.rightDensity;
        points[i].velocity_tau = borderSoda.rightVelocity;
        points[i].velocity_normal = 0;
        points[i].velocity = borderSoda.rightVelocity;
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
    }
    system->prepareSolving(points);
}

void GodunovSolverSoda::updatePoints()
{
    auto size = points.size();
    for(size_t i = 0; i < size; i++)
    {
        points[i].velocity_tau = system->getVelocityTau(i);
        points[i].velocity_normal = system->getVelocityNormal(i);
        points[i].velocity = system->getVelocity(i);
        points[i].density = system->getDensity(i);
        points[i].pressure = system->getPressure(i);
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
    }
}



