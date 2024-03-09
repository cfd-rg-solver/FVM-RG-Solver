//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
#include <filesystem>


std::string GetCurrentWorkingDir( void ) {
    std::filesystem::path currentWorkingDir = std::filesystem::current_path();
    std::filesystem::path parentDir = currentWorkingDir.parent_path().parent_path();

   std::string res = parentDir.string() + "/FVM-RG-Solver/example-couette";

    // std::string res = parentDir.string() + "/main/example-couette"; // !normal case using qt
    return res;
//    return currentWorkingDir; // without qt
//    return currentWorkingDir.string();
}

namespace fs = std::filesystem;
int main()
{

    std::string outputData = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << outputData << std::endl;

    //////////////////////////////////////////////////////////////
    ///////////////////// Border Condition for Couette ///////////
    //////////////////////////////////////////////////////////////
    double T_up_wall = 1000;
    double T_down_wall = 1000;
    double velocity_up = 300;
    double velocity_down = 0;

    BorderConditionCouette borderConditionCouette;
    borderConditionCouette.setWallParameters(velocity_up, velocity_down, T_up_wall, T_down_wall);

    //////////////////////////////////////////////////////////////
    ///////////////////// Border Condition for Couette ///////////
    /////////////////////////// Personal /////////////////////////
    T_up_wall = 1000;
    velocity_up = 0;
    double start_velocity_normal = 10;

    BorderConditionPersonal borderConditionPersonal;
    borderConditionPersonal.setWallParameters(0, 0, T_up_wall, T_down_wall); // в данном кейсе (со сплошной среджой и стенкой) значимым является только значения на верхней стенке

    //////////////////////////////////////////////////////////////


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
    ////////////////////////////  O2_O  /////////////////////////

    UniformDistributionBorder startParamCouetteO2_O;
    UniformDistributionBorderPersonal startParamCouetteO2_OPersonal;
    macroParam startParamO2_O(O2_O);
    startParamO2_O.density = 0.03168;
    startParamO2_O.fractionArray[0] = 0.99;
    startParamO2_O.densityArray[0] =  startParamO2_O.fractionArray[0] * startParamO2_O.density;

    startParamO2_O.fractionArray[1] = 0.01;
    startParamO2_O.densityArray[1] =  startParamO2_O.fractionArray[1] * startParamO2_O.density;

    startParamO2_O.temp = 140; //140
    startParamO2_O.velocity_tau = 0;
    startParamO2_O.velocity_normal = 0;

    startParamCouetteO2_O.setBorderCondition(&borderConditionCouette);
    startParamCouetteO2_O.setDistributionParameter(startParamO2_O);

    startParamCouetteO2_OPersonal.setBorderCondition(&borderConditionPersonal);
    startParamCouetteO2_OPersonal.setDistributionParameter(startParamO2_O);
    startParamCouetteO2_OPersonal.setNormalVelocity(start_velocity_normal);
    //////////////////////////////////////////////////////////////
    ///////////////////// Start param for Couette ////////////////
    ////////////////////////////  Ar  ///////////////////////////

    UniformDistributionBorder startParamCouetteAr;
    UniformDistributionBorderPersonal startParamCouetteArPersonal;
    macroParam startParamAr(Ar);
    startParamAr.density = 0.03168;
    startParamAr.fractionArray[0] = 1;
    startParamAr.densityArray[0] =  startParamAr.fractionArray[0] * startParamAr.density;

    startParamAr.temp = 140; //140
    startParamAr.velocity_tau = 0;
    startParamAr.velocity_normal = 0;


    startParamCouetteAr.setBorderCondition(&borderConditionCouette);
    startParamCouetteAr.setDistributionParameter(startParamAr);

    startParamCouetteArPersonal.setBorderCondition(&borderConditionPersonal);
    startParamCouetteArPersonal.setDistributionParameter(startParamAr);
    startParamCouetteArPersonal.setNormalVelocity(start_velocity_normal);


    //////////////////////////////////////////////////////////////

    solverParams solParam;
    solParam.NumCell     = 202;    // Число расчтеных ячеек с учетом двух фиктивных ячеек
//    solParam.Gamma    = 1.67;   // Ar
    solParam.Gamma    = 1.32;   // O2_O
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

//    GodunovSolver solver(Ar ,solParam, SystemOfEquationType::couette2Alt, RiemannSolverType::HLLESolver);
    GodunovSolver solver(O2_O ,solParam, SystemOfEquationType::couette2AltBinary, RiemannSolverType::HLLESolver);
    double h = 1;
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));


    solver.setBorderConditions(&borderConditionCouette);  // for couette

    solver.setStartDistribution(&startParamCouetteO2_O); // for couette Ar

    solver.solve();
}
