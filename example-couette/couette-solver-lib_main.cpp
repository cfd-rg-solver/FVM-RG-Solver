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

//#define Liia
namespace fs = std::filesystem;
int main()
{

    std::string outputData = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << outputData << std::endl;

    MixtureComponent argon;
    argon.name = "Ar";
    argon.density = 1.7839e-3; // case 1
    //argon.density = 0.800773;
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21/kB; //! Mistake
    argon.sigma = 3.33E-10;
    Mixture mixture({argon});

    macroParam startParam(mixture);
    startParam.fractionArray[0] = 1;
    //startParam.pressure = 218563.81; //218563.81
    startParam.densityArray[0] = argon.density;
    startParam.density = argon.density;
    startParam.temp = 1000; //140
    startParam.velocity_tau = 0.00001;

    solverParams solParam;
    solParam.NumCell     = 102;    // Число расчтеных ячеек с учетом двух фиктивных ячеек
    solParam.Gamma    = 1.67;    // Показатель адиабаты
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
    GodunovSolver solver(mixture,startParam,solParam, SystemOfEquationType::couette2Alt, RiemannSolverType::HLLESolver);
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));
    solver.setBorderConditions(velocity,T2wall,T1wall);
    solver.solve();
}
