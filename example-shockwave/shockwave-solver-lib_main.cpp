//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
#include <filesystem>


std::string GetCurrentWorkingDir( void ) {
    std::filesystem::path currentWorkingDir = std::filesystem::current_path();
    std::filesystem::path parentDir = currentWorkingDir.parent_path().parent_path();

    std::string res = parentDir.string() + "/FVM-RG-Solver/example-shockwave"; // ! win case

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

    // swtting Ar mixture
    MixtureComponent argon;
    argon.name = "Ar";
    //argon.density = 0.03168*0.99; // initial concentration - 0.99, pressure 133 000 Pa
    //argon.density = 0.800773;
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21/kB; //! Mistake
    argon.sigma = 3.33E-10;
    argon.numberAtoms = 1;
    argon.omega_e = 0;

    std::vector<MixtureComponent> tmp2 = {argon};
    Mixture Ar(tmp2);

    //////////////////////////////////////////////////////////////
    ///////////////////// Border Condition for Couette ///////////
    //////////////////////////////////////////////////////////////
    double T_up_wall = 1000;
    double T_down_wall = 1000;
    double velocity_up = 300;
    double velocity_down = 0;

    BorderConditionCouette borderConditionCouette;
    borderConditionCouette.setWallParameters(velocity_up, velocity_down, T_up_wall, T_down_wall);


    ///////////////////////////////////////////////////////////////
    ///////////////////// Border Condition for Shock Wave ////////
    /////////////////////////////////////////////////////////////
    ///
    // рассматриваем уравнения граничных условий, пусть left = 0, right = n:
    double velocity_left = 1227.4; // shock wave, so the velocity is supersonic, let's set it to 1615 m/s ~ 5 Ma for argon at room temperature
    double density_left = 0.00010666; // kg/m^3, calculated for atmospheric pressure
    double T_left = 300; // Kelvin
    double pressure_left = UniversalGasConstant * T_left * density_left / argon.molarMass;

    double velocity_right = 370.7071;
    double density_right = 0.000353;
    double T_right = 1612.93;
    double pressure_right = UniversalGasConstant * T_right * density_right / argon.molarMass;

    BorderConditionShockwave borderConditionShockwave;
    borderConditionShockwave.setBorderParameters(
        velocity_left, density_left, T_left,
        velocity_right, density_right, T_right
        );


    //////////////////////////////////////////////////////////////
    ///////////////////// Start param for Couette ////////////////
    ////////////////////////////  Ar  ///////////////////////////

    UniformDistributionBorder startParamCouetteAr;
    // UniformDistributionBorderPersonal startParamCouetteArPersonal;
    macroParam startParamAr(Ar);
    startParamAr.density = 0.03168;
    startParamAr.fractionArray[0] = 1;
    startParamAr.densityArray[0] =  startParamAr.fractionArray[0] * startParamAr.density;

    startParamAr.temp = 140; //140
    startParamAr.velocity_tau = 0;
    startParamAr.velocity_normal = 0;


    startParamCouetteAr.setBorderCondition(&borderConditionCouette);
    startParamCouetteAr.setDistributionParameter(startParamAr);

    //////////////////////////////////////////////////////////////
    ////////////////// Start param for Shockwave /////////////////
    ////////////////////////////  Ar  ///////////////////////////
    GapDistribution startParamShockwaveAr;
    macroParam leftStartParam(Ar);
    macroParam rightStartParam(Ar);

    leftStartParam.density = density_left;
    leftStartParam.pressure = pressure_left;
    leftStartParam.velocity = velocity_left;
    leftStartParam.fractionArray[0] = 1;
    leftStartParam.densityArray[0] =  leftStartParam.fractionArray[0] * leftStartParam.density;

    rightStartParam.density = density_right;
    rightStartParam.pressure = pressure_right;
    rightStartParam.velocity = velocity_right;
    rightStartParam.fractionArray[0] = 1;
    rightStartParam.densityArray[0] = rightStartParam.fractionArray[0] * rightStartParam.density;

    startParamShockwaveAr.setDistributionParameter(leftStartParam, rightStartParam);

    //////////////////////////////////////////////////////////////
    ///
    ///

    double viscocity_argon = 22.7e-5; // approximate argon viscocity at low prassure
    double MFP = viscocity_argon / pressure_left * sqrt(M_PI * UniversalGasConstant * T_left / argon.molarMass); // mean free path length

    std::cout << "mean free path: " << MFP << std::endl;

    solverParams solParam;
    solParam.NumCell     = 1e2 + 2;    // Число расчтеных ячеек с учетом двух фиктивных ячеек
    solParam.Gamma    = 1.67;   // Ar
    // solParam.Gamma    = 1.32;   // O2_O
    solParam.CFL      = 0.9;    // Число Куранта
    solParam.MaxIter     = 1000000; // максимальное кол-во итареций
    solParam.Ma       = 3.8;    // Число Маха

    double precision = 1E-7; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(10000);


    DataWriter writer(outputData);
    DataReader reader(outputData + "\prev_data");

    reader.read();
    vector<macroParam> startParameters;
    reader.getPoints(startParameters);

    // GodunovSolver solver(Ar ,solParam, SystemOfEquationType::couette2Alt, RiemannSolverType::HLLESolver);
    GodunovSolver solver(Ar ,solParam, SystemOfEquationType::shockwave1, RiemannSolverType::HLLESolver);
    double h = 30 * MFP; // m
    // double h = 1;
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));


    // solver.setBorderConditions(&borderConditionCouette);  // for couette
    solver.setBorderConditions(&borderConditionShockwave);

    // solver.setStartDistribution(&startParamCouetteAr); // for couette Ar
    solver.setStartDistribution(&startParamShockwaveAr);

    solver.solve();
}
