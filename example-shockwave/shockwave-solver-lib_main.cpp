#include "godunovsolver.h"
#include "bordercondition.h"
#include "mixture.h"
#include <filesystem>
#include <iostream>


std::string GetCurrentWorkingDir( void ) {
    std::filesystem::path currentWorkingDir = std::filesystem::current_path();
    std::filesystem::path parentDir = currentWorkingDir.parent_path().parent_path();
    std::string res = parentDir.string() + "/FVM-RG-Solver/example-shockwave";
    return res;
}

namespace fs = std::filesystem;
int main()
{
    std::string outputData = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << outputData << std::endl;

    //////////////////////////////////////////////////////////////
    ///////////////// Border Condition for Shockwave ////////////
    //////////////////////////// CH4 ////////////////////////////

    /// CH4 (methane) /// 
    //MixtureComponent methane;
    //methane.name = "methane";
    //methane.molarMass = 0.0160425; // kg/mol
    //methane.mass = 2.656e-26; // kg
    //methane.numberAtoms = 5;
    // methane.density = 0.657; // kg/m^3, when about 20 C
    // methane.numberVibrLvl = ?;
    // methane.epsilonDevK = ?
    // methane.sigma = ?
    // methane.omega_e = ?
    /////////////////////

    /// Create mixture ///
    // std::vector<MixtureComponent> tmp1 = { methane }; // - однокомпонентная смесь
    // Mixture CH4(tmp1);
    /////////////////////

    //////////////////////////////////////////////////////////////
    ////////////////// Start param for Shockwave /////////////////
    ////////////////////////////  CH4  ///////////////////////////
    // UniformDistributionBorder startParamShockwaveCH4;
    // UniformDistributionBorderPersonal startParamShockwaveCH4Personal;
    // macroParam startParamCH4(CH4);

    // startParamCH4.temp = 200; // Kelvin
    // startParamCH4.velocity_normal = 0;

    // startParamCH4.density = 0.971; // kg/m^3, 200K https://www.engineeringtoolbox.com/methane-density-specific-weight-temperature-pressure-d_2020.html
    // startParamCH4.fractionArray[0] = 1; // однокомпонентная смесь
    // startParamCH4.densityArray[0] = startParamCH4.fractionArray[0] * startParamCH4.density;

    // startParamShockwaveCH4.setBorderCondition(&borderConditionShockwave);
    // startParamShockwaveCH4.setDistributionParameter(startParamCH4);
    //////////////////////////////////////////////////////////////

    
    //////////////////////////////////////////////////////////////
    ///////////////// Border Condition for Shockwave ////////////
    //////////////////////////// Ar ////////////////////////////

    /// Ar (argon) /// 
    MixtureComponent argon;
    argon.name = "argon";
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992e-26;
    argon.numberAtoms = 1;
    // argon.numberVibrLvl = ?;
    // argon.epsilonDevK = ?
    // argon.sigma = 3.33e-10; ?
    // argon.omega_e = ?
    /////////////////////

    // рассматриваем уравнения граничных условий,
    // пусть left = 0, right = n:
    double velocity_left = 1615; // shock wave, so the velocity is supersonic, let's set it to 1615 m/s ~ 5 Ma for argon at room temperature
    double density_left = 1.759942; // kg/m^3, calculated for atmospheric pressure
    double T_left = 300; // Kelvin
    double pressure_left = UniversalGasConstant * T_left * density_left / argon.molarMass;
    // double energy_left = 3 * kB * T_left / (2 * argon.mass); // одноатомный аргон

    double velocity_right = 452.0778769113248; // recalculated ...
    double density_right = 6.287205092669588;
    double T_right = 2609.9281339791883; // smth is wrong, even for initial conditions it is a bit strange that the gas is so cold on the right border
    double pressure_right = UniversalGasConstant * T_right * density_right / argon.molarMass;

    BorderConditionShockwave borderConditionShockwave;
    borderConditionShockwave.setBorderParameters(
        velocity_left, density_left, T_left,
        velocity_right, density_right, T_right);
    /////////////////////

    /// Create mixture ///
    std::vector<MixtureComponent> tmp2 = { argon }; // - однокомпонентная смесь
    Mixture Ar(tmp2);
    /////////////////////

    //////////////////////////////////////////////////////////////
    ////////////////// Start param for Shockwave /////////////////
    ////////////////////////////  Ar  ///////////////////////////
    GapDistribution startParamShockwaveAr;
    macroParam leftStartParam(Ar);
    macroParam rightStartParam(Ar);
    leftStartParam.density = density_left;
    leftStartParam.pressure = pressure_left;
    leftStartParam.velocity = velocity_left;
    rightStartParam.density = density_right;
    rightStartParam.pressure = pressure_right;
    rightStartParam.velocity = velocity_right;
    startParamShockwaveAr.setDistributionParameter(leftStartParam, rightStartParam);
    //////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////
    /////////////////////// General Solver ///////////////////////
    //////////////////////////////////////////////////////////////
    solverParams solParam;
    solParam.NumCell = 102; // число ячеек (с учетом двух фиктивных)
    solParam.Gamma = 1.667; // показатель адиабаты c_p/c_v, argon
    solParam.CFL = 0.9; // число Куранта
    solParam.MaxIter = 1e7; // максимальное кол-во итераций
    solParam.Ma = 5; // число Маха

    double precision = 1e-5; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(1e4);
    //////////////////////////////////////////////////////////////   

    ////////////////////////////////////////////////////////////// 
    DataWriter writer(outputData);
    DataReader reader(outputData + "\prev_data");

    reader.read();
    vector<macroParam> startParameters;
    reader.getPoints(startParameters);
    //////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////
    // GodunovSolver solver(CH4, solParam, SystemOfEquationType::shockwave2, RiemannSolverType::HLLESolver);
    GodunovSolver solver(Ar, solParam, SystemOfEquationType::shockwave1, RiemannSolverType::HLLESolver);
    double h = 1;
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));
    solver.setBorderConditions(&borderConditionShockwave);  // for shockwave
    solver.setStartDistribution(&startParamShockwaveAr); // for shockwave
    solver.solve();
    //////////////////////////////////////////////////////////////

}
