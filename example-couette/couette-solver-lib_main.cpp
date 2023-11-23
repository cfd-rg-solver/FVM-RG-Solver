//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
#include <filesystem>


std::string GetCurrentWorkingDir( void ) {
    std::filesystem::path currentWorkingDir = std::filesystem::current_path();
    std::filesystem::path parentDir = currentWorkingDir.parent_path().parent_path();

    std::string res = parentDir.string() + "/FVM-RG-Solver/main/example-couette";

//    std::string res = parentDir.string() + "/main/example-couette"; // !normal case using qt
    return res;
//    return currentWorkingDir; // without qt
//    return currentWorkingDir.string();
}

namespace fs = std::filesystem;
int main()
{

    std::string outputData = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << outputData << std::endl;

    // Ar

    MixtureComponent argon;
    argon.name = "Ar";
    argon.density = 0.03168*0.99; // initial concentration - 0.99, pressure 133 000 Pa
    //argon.density = 0.800773;
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21/kB; //! Mistake
    argon.sigma = 3.33E-10;
    Mixture Ar({argon});

    // O2/O

    MixtureComponent Oxygen2;
    Oxygen2.name = "O2";
//    Oxygen2.density = 1.42897e-3*0.99;
    Oxygen2.density = 0.03168;
    Oxygen2.molarMass = 0.0319988;
    Oxygen2.mass = 5.313725014956E-26;
    Oxygen2.epsilonDevK = 1.48281651048e-21/kB;
    Oxygen2.sigma = 3.5155e-10;
    Oxygen2.numberAtoms = 2;
    Oxygen2.numberVibrLvl = 15;
    Oxygen2.omega_e = 1580.19; // см^-1

    MixtureComponent Oxygen;
    Oxygen.name = "O";
    Oxygen.density = 0.03168; // initial concentration - 0.01
    Oxygen.molarMass = 0.015999;
    Oxygen.mass = 2.6567628316576e-26;
    Oxygen.epsilonDevK = 1.1045188159999998e-21/kB;
    Oxygen.sigma = 2.75E-10;
    Oxygen.numberAtoms = 1;

    std::vector<MixtureComponent> tmp = {Oxygen2,Oxygen};
    Mixture O2_O(tmp);

    /////Выбор смеси/////

    Mixture mixture = O2_O;

    /////////////////////

    macroParam startParam(mixture);
    startParam.fractionArray[0] = 0.99;
    startParam.densityArray[0] =  startParam.fractionArray[0] * mixture.components[0].density;

    startParam.fractionArray[1] = 0.01;
    startParam.densityArray[1] =  startParam.fractionArray[1] * mixture.components[1].density;

    startParam.density = 0;
    for(int i = 0; i < mixture.NumberOfComponents; i++)
        startParam.density += startParam.densityArray[i];
    startParam.temp = 800; //140
    startParam.velocity_tau = 0.00001;


    solverParams solParam;
    solParam.NumCell     = 102;    // Число расчтеных ячеек с учетом двух фиктивных ячеек
    solParam.Gamma    = 1.32;
//    solParam.Gamma    = 1.67;    // Показатель адиабаты, Ar
    solParam.CFL      = 0.9;    // Число Куранта
    solParam.MaxIter     = 10000000; // максимальное кол-во итареций
    solParam.Ma       = 0.1;    // Число маха

    double precision = 1E-5; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(10000);


    DataWriter writer(outputData);
    DataReader reader(outputData + "\prev_data");

    reader.read();
    vector<macroParam> startParameters;
    reader.getPoints(startParameters);
    double T1wall = 1000;
    double T2wall = 1000;
    double velocity = 300;
    double h = 1;
    GodunovSolver solver(mixture,startParam,solParam, SystemOfEquationType::couette2AltBinary, RiemannSolverType::HLLESolver);
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));
    solver.setBorderConditions(velocity,T2wall,T1wall);
    solver.solve();
}
