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

        system->computeF(points, delta_h);
        if(system->systemType == SystemOfEquationType::couette2Alt || system->systemType == SystemOfEquationType::couette2AltBinary)
            system->computeFv(points, delta_h);
        riemannSolver->computeFlux(system);
        //riemannSolver->computeFlux(system, delta_h);
        //riemannSolver->computeFlux(system,timeSolvind.last(),delta_h);

        // Вычисляем вектор релаксационных членов
        //computeR();

        // Обновляем вектор U
        system->updateU(delta_h,timeSolvind.last());

        cout<<"iter ---------  "<< i<<endl;
        cout<<points[0].density<<" "<<points[0].temp<<" "<<points[0].pressure<<" "<<points[0].densityArray[0]<<endl;
        cout<<system->F[0][0]<<" "<<system->F[1][0]<<" "<<system->F[2][0]<<" "<<system->F[3][0]<<" "<<system->F[4][0]<<endl;
        cout<<system->F[0][1]<<" "<<system->F[1][1]<<" "<<system->F[2][1]<<" "<<system->F[3][1]<<" "<<system->F[4][1]<<endl;
        cout<<system->Flux[0][1]<<" "<<system->Flux[1][1]<<" "<<system->Flux[2][1]<<" "<<system->Flux[3][1]<<" "<<system->Flux[4][1]<<endl;
        cout<<system->U[0][0]<<" "<<system->U[1][0]<<" "<<system->U[2][0]<<" "<<system->U[3][0]<<" "<<system->U[4][0]<<endl;
        cout<<system->U[0][1]<<" "<<system->U[1][1]<<" "<<system->U[2][1]<<" "<<system->U[3][1]<<" "<<system->U[4][1]<<endl<<endl;

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
    }
    writePoints(T*1000000);
}
