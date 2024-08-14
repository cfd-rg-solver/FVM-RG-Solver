#include "godunovsolver.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <omp.h>

void GodunovSolver::setOutputDirectory(const std::string& outputDir) {
    // ! later connect it with writer
    outputDirectory = outputDir;
    // Ensure the directory exists
    std::filesystem::create_directory(outputDirectory);
}

void GodunovSolver::solve()
{
    lawChecker.setDh(delta_h);
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
            // std::cout << "energy " << system->getEnergy(30) << std::endl;
            std::cout << "max wave speed " << max << std::endl;
        }
        //проверка точности
        if(isObserverWatching)
        {
            // то есть если проверка наблюдателя не пройдена, нужно прекратить рассчёт
            if(!observerCheck(i))
                break;
        }
        lawCheck(i);
    }
    writePoints(T*1000000);
    
    // Print viscous flux components
    printViscousFluxComponents(system, points.size());
}


void GodunovSolver::printViscousFluxComponents(SystemOfEquation* system, int N) {
    // Open files for writing
    std::ofstream stressTensorFile(outputDirectory + "/data/stress_tensor.txt");
    std::ofstream heatFluxFile(outputDirectory + "/data/heat_flux.txt");

    if (!stressTensorFile.is_open() || !heatFluxFile.is_open()) {
        std::cerr << "Error opening files for writing." << std::endl;
        return;
    }

    stressTensorFile << 'y' << " "  <<  'P' << std::endl;
    heatFluxFile << 'y' << " "  << 'q' << std::endl;

    // Iterate over the calculation area
    for (int i = 0; i < N; i++) {
        // Extract the stress tensor component
        // std::cout << "i: " << i << std::endl;
        double P_i = system->getNormalStress(i);
        // std::cout << "P_i: " << P_i << std::endl;
        // Extract the heat flux component
        double q_i = system->getHeatFlux(i);
        // std::cout << "q_i: " << q_i << std::endl;
        
        // Write to the respective files
        stressTensorFile << delta_h * i  << " "  << P_i << std::endl;
        heatFluxFile << delta_h * i  << " " << q_i << std::endl;
    }

    // Close the files
    stressTensorFile.close();
    heatFluxFile.close();
}
