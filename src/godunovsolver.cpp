#include "godunovsolver.h"
#include "lawchecker.h"

#include <iostream>
#include <omp.h>
void GodunovSolver::solve()
{
    LawChecker cheker(0,201,100000,100001,delta_h);
    writePoints(-1);
    double T = 0;
    for(size_t i  = 0; i < solParam.MaxIter; i++)
    {
        // Устанавливаем текущий временной шаг
        setDt();
        T += timeSolvind.last();

        system->computeF(points, delta_h);
        if(system->systemType == SystemOfEquationType::couette2Alt || system->systemType == SystemOfEquationType::couette2AltBinary)
            system->computeFv(points, delta_h);
        //riemannSolver->computeFlux(system);
        //riemannSolver->computeFlux(system, delta_h);
        riemannSolver->computeFlux(system,timeSolvind.last(),delta_h);

        // Вычисляем вектор релаксационных членов
        //computeR();


        // Обновляем вектор U
        system->updateU(delta_h,timeSolvind.last());

        // Обновляем вектор макропараметров
        updatePoints();
        // обновляем вектор U с учётом граничных условий
        system->updateBorderU(points);

        //записать данные, если это требуется
        //writePoints(T*1000000); // микросек

        double max;
        if(i%10000 == 0)
        {
            std::cout<<i<<" iteration"<<std::endl;
            writePoints(T*1000000); // микросек

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
        cheker.rememberU(system->U[system->energy],i,timeSolvind.last());
        cheker.rememberF(system->F[system->energy], system->Fv[system->energy],i,timeSolvind.last());
    }
    writePoints(T*1000000);
}
