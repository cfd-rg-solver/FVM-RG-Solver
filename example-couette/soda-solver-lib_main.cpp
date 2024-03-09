//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "datawriter.h"
#include "observer.h"
#include <filesystem>

namespace fs = std::filesystem;
int main()
{
    MixtureComponent argon;
    argon.name = "Ar";
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21;
    argon.sigma = 3.33E-10;
    Mixture Ar({argon});

    macroParam startParam(Ar);
    startParam.density = 1.7839;
    startParam.fractionArray[0] = 1;
    startParam.pressure = 218563.81; //218563.81
    //startParam.densityArray[0] = argon.density;
    startParam.temp = 273; //140

    solverParams solParam;
    solParam.NumCell     = 100;    // Число расчтеных ячеек
    solParam.Gamma    = 1.4;    // Показатель адиабаты
    solParam.CFL      = 0.8;    // Число Куранта
    solParam.MaxIter     = 1000; // максимальное кол-во шагов по времени
    solParam.Ma       = 0.5;    // Число маха

    double precision = 0.000001; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(1000);

    // это меняешь под себя. Он так создаст папку data
    // если не использовать setWriter, то записи не будет, но папка создастся, ибо она в конструкторе зашита
    // он автоматически очищает папку перед новым рассчётом
    DataWriter writer("D:/couette/couette-solver-lib");

    BorderConditionSoda borderSoda;
    borderSoda.setBorderParameters( 0,  1,  1,
                             0,  0.125,  0.1);

    GapDistribution distributionSoda;
    macroParam leftStartParam(Ar);
    macroParam rightStartParam(Ar);

    leftStartParam.density = 1;
    leftStartParam.pressure = 1;
    leftStartParam.velocity_tau = 0;
    leftStartParam.fractionArray[0] = 1;
    leftStartParam.densityArray[0] =  leftStartParam.fractionArray[0] * leftStartParam.density;

    rightStartParam.density = 0.125;
    rightStartParam.pressure = 0.1;
    rightStartParam.velocity_tau = 0;
    rightStartParam.fractionArray[0] = 1;
    rightStartParam.densityArray[0] = rightStartParam.fractionArray[0] * rightStartParam.density;

    distributionSoda.setDistributionParameter(leftStartParam, rightStartParam);

    double h = 1;

    GodunovSolver solver(Ar,solParam, SystemOfEquationType::soda,RiemannSolverType::HLLESolver);
    solver.setDelta_h(h / solParam.NumCell);
    solver.setWriter(&writer);
    //solver.setObserver(&watcher);
    solver.setBorderConditions(&borderSoda);
    solver.setStartDistribution(&distributionSoda);
    solver.solve();
}
