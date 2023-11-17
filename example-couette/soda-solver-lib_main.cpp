//#include "hllcsolver.h"
#include "godunovsolversoda.h"
#include "datawriter.h"
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
    argon.density = 1.7839;
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21;
    argon.sigma = 3.33E-10;
    Mixture mixture({argon});

    macroParam startParam(mixture);
    startParam.fractionArray[0] = 1;
    startParam.pressure = 218563.81; //218563.81
    //startParam.densityArray[0] = argon.density;
    startParam.temp = 273; //140

    solverParams solParam;
    solParam.NumCell     = 100;    // Число расчтеных ячеек
    solParam.Gamma    = 1.4;    // Показатель адиабаты
    solParam.CFL      = 0.8;    // Число Куранта
    solParam.MaxIter     = 100000; // максимальное кол-во шагов по времени
    solParam.Ma       = 0.5;    // Число маха

    double precision = 0.000001; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(1000);

    // это меняешь под себя. Он так создаст папку data
    // если не использовать setWriter, то записи не будет, но папка создастся, ибо она в конструкторе зашита
    // он автоматически очищает папку перед новым рассчётом
    DataWriter writer(outputData);
    DataReader reader(outputData + "\prev_data");

    BorderConditionSoda borderSoda;
    borderSoda.leftDensity = 1;
    borderSoda.leftPressure = 1;
    borderSoda.leftVelocity = 0;
    borderSoda.rightDensity = 0.125;
    borderSoda.rightPressure = 0.1;
    borderSoda.rightVelocity = 0;

    double h = 1;

    GodunovSolverSoda solver(mixture,startParam,solParam, SystemOfEquationType::soda,RiemannSolverType::HLLCSolver);
    solver.setDelta_h(h / solParam.NumCell);
    solver.setWriter(&writer);
    //solver.setObserver(&watcher);
    solver.setBorderConditions(borderSoda);
    solver.solve();
}
