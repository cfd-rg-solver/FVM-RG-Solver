#include "godunovsolver.h"
#include <iostream>
#include <time.h>
#include <omp.h>
void GodunovSolver::solve()
{
    writePoints(-1);
    double T = 0;
    clock_t start = clock();
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
            system->systemType == SystemOfEquationType::shockwave1 ||
            system->systemType == SystemOfEquationType::shockwave2)
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
        // writePoints(T*1000000); // микросек



        double max;
        if(i%100 == 0)
        {
            std::cout<<i<<" iteration"<<std::endl;
            writePoints(T*1000000); // микросек

            max = riemannSolver->maxSignalVelocity;
            std::cout << "energy " << system->getEnergy(30) << std::endl;
            std::cout << "max wave speed " << max << std::endl;
            clock_t end = clock();
            double seconds = (double)(end - start) / CLOCKS_PER_SEC;
            printf("Time of solution till current iter: %f seconds\n", seconds);
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
