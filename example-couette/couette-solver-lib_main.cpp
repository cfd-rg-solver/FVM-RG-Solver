//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
#include <filesystem>


std::string GetCurrentWorkingDir( void ) {
    std::filesystem::path currentWorkingDir = std::filesystem::current_path();
    std::filesystem::path parentDir = currentWorkingDir.parent_path().parent_path();

//    std::string res = parentDir.string() + "/FVM-RG-Solver/main/example-couette"; ! win case

    std::string res = parentDir.string() + "/main/example-couette"; // !normal case using qt
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
    //argon.density = 0.03168*0.99; // initial concentration - 0.99, pressure 133 000 Pa
    //argon.density = 0.800773;
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21/kB; //! Mistake
    argon.sigma = 3.33E-10;

    // O2/O

    MixtureComponent Oxygen2;
    Oxygen2.name = "O2";
//    Oxygen2.density = 1.42897e-3*0.99;
    Oxygen2.molarMass = 0.0319988;
    Oxygen2.mass = 5.313725014956E-26;
    Oxygen2.epsilonDevK = 1.48281651048e-21/kB;
    Oxygen2.sigma = 3.5155e-10;
    Oxygen2.numberAtoms = 2;
    Oxygen2.numberVibrLvl = 15;
    Oxygen2.omega_e = 1580.19; // см^-1

    MixtureComponent Oxygen;
    Oxygen.name = "O";
    Oxygen.molarMass = 0.015999;
    Oxygen.mass = 2.6567628316576e-26;
    Oxygen.epsilonDevK = 1.1045188159999998e-21/kB;
    Oxygen.sigma = 2.75E-10;
    Oxygen.numberAtoms = 1;

    std::vector<MixtureComponent> tmp1 = {Oxygen2,Oxygen};
    Mixture O2_O(tmp1);

    std::vector<MixtureComponent> tmp2 = {argon};
    Mixture Ar(tmp2);
    //////////////////////////////////////////////////////////////
    ///////////////////// Start param for Couette ////////////////
    //////////////////////////////////////////////////////////////

    macroParam startParam(O2_O);
    startParam.density = 0.03168;
    startParam.fractionArray[0] = 0.99;
    startParam.densityArray[0] =  startParam.fractionArray[0] * startParam.density;

    startParam.fractionArray[1] = 0.01;
    startParam.densityArray[1] =  startParam.fractionArray[1] * startParam.density;

    startParam.temp = 800; //140
    startParam.velocity_tau = 0.00001;

    //////////////////////////////////////////////////////////////
    ///////////////////// Start param for Soda ///////////////////
    //////////////////////////////////////////////////////////////

    macroParam leftStartParam(Ar);
    macroParam rightStartParam(Ar);

    leftStartParam.density = 1;
    leftStartParam.pressure = 1;
    leftStartParam.velocity = 0;
    rightStartParam.density = 0.125;
    rightStartParam.pressure = 0.1;
    rightStartParam.velocity = 0;

    //////////////////////////////////////////////////////////////

    solverParams solParam;
    solParam.NumCell     = 202;    // Число расчтеных ячеек с учетом двух фиктивных ячеек
    solParam.Gamma    = 1.67;
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
    GodunovSolver solver(O2_O,solParam, SystemOfEquationType::couette2AltBinary, RiemannSolverType::HLLESolver);
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));

    // в следующих строчках супер сложные зависимости, скорее всего надо создавать какой-то отдельный класс для начальных и граничных уловий,
    // ибо то, что есть сейчас супер завязано на определнной последовательности и выглядит не интуитивно
    // т.е. что-то задается в конструкторе, что-то в вызыываемом методе, что-то не может быть вызвано без предварительного вызова другого и т.д, хотя логически
    // оно немного и связано (типо начальное распределение, а точнее его значение в фиктивных ячейках, зависит от граничных условий),
    solver.setBorderConditions(velocity,T2wall,T1wall);

    solver.setStartDistribution(startParam); // for couette
    //solver.setStartDistribution(leftStartParam, rightStartParam); // for soda

    solver.solve();
}
