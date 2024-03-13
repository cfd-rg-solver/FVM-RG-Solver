#include "godunovsolver.h"
#include <iostream>
#include <omp.h>
void GodunovSolver::solve()
{
    writePoints(-1);
    double T = 0;
    for(size_t i  = 0; i < solParam.MaxIter; i++)
    {
        // Устанавливаем текущий временной шаг
        setDt();
        T += timeSolvind.last();
        // if(i%10000 == 0)
        // {
        //     std::cout<<i<<" next time set "<< T <<std::endl;
        // }

        system->computeF(points, delta_h);

        if(system->systemType == SystemOfEquationType::couette2Alt ||
            system->systemType == SystemOfEquationType::couette2AltBinary ||
            system->systemType == SystemOfEquationType::shockwave1)
            system->computeFv(points, delta_h);

        riemannSolver->computeFlux(system);

        //riemannSolver->computeFlux(system, delta_h);
        //riemannSolver->computeFlux(system,timeSolvind.last(),delta_h);

        // Вычисляем вектор релаксационных членов
        //computeR();

        // Обновляем вектор U
        system->updateU(delta_h,timeSolvind.last());

        // Обновляем вектор макропараметров
        updatePoints();

        // обновляем вектор U с учётом граничных условий
        system->updateBorderU(points); // this one should be calculated on the basis of boundary conditions type

        //записать данные, если это требуется
        //writePoints(T*1000000); // микросек

        // std::cout << "visc " << coeffSolver->shareViscositySimple(points[(int)(solParam.NumCell/2)]) << std::endl;

        double max;
        if(i%1000 == 0)
        {
            std::cout<<i<<" iteration"<<std::endl;
            writePoints(i); // микросек

            max = riemannSolver->maxSignalVelocity;
            std::cout << "max wave speed " << max << std::endl;

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
