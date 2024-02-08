#include "hdr/godunovsolver.h"
#include "hdr/bordercondition.h"
#include <filesystem>


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
    /////////////////// Border Condition for Shockwave ///////////
    //////////////////////////////////////////////////////////////
    double T_up_wall = 1000;
    double velocity_up = 0;
    double pressure_up = 1; // ??? todo 

    BorderConditionShockwave borderConditionShockwave;
    borderConditionShockwave.setWallParameters(velocity_up, T_up_wall, pressure_up);

    /// CH4 (methane) /// 
    MixtureComponent methane;
    methane.name = "methane";
    methane.molarMass = 0.0160425; // kg/mol
    methane.mass = 2.656e-26; // kg
    methane.numberAtoms = 5;
    // methane.density = 0.657; // kg/m^3, when about 20 C
    // methane.numberVibrLvl = ?;
    // methane.epsilonDevK = ?
    // methane.sigma = ?
    // methane.omega_e = ?
    /////////////////////

    /// Create mixture ///
    std::vector<MixtureComponent> tmp = {methane}; // - однокомпонентная смесь
    Mixture CH4(tmp);
    /////////////////////

    //////////////////////////////////////////////////////////////
    ////////////////// Start param for Shockwave /////////////////
    ////////////////////////////  CH4  ///////////////////////////
    UniformDistributionBorder startParamShockwaveCH4;
    UniformDistributionBorderPersonal startParamShockwaveCH4Personal;
    macroParam startParamCH4(CH4);

    startParamCH4.temp = 200; // Kelvin
    startParamCH4.velocity_tau = 0;
    startParamCH4.velocity_normal = 0;

    startParamCH4.density = 0.971; // kg/m^3, 200K https://www.engineeringtoolbox.com/methane-density-specific-weight-temperature-pressure-d_2020.html
    startParamCH4.fractionArray[0] = 1; // однокомпонентная смесь
    startParamCH4.densityArray[0] =  startParamCH4.fractionArray[0] * startParamCH4.density;

    startParamShockwaveCH4.setBorderCondition(&borderConditionShockwave);
    startParamShockwaveCH4.setDistributionParameter(startParamCH4);
    //////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////
    solverParams solParam;
    solParam.NumCell = 102; // число ячеек (с учетом двух фиктивных)
    solParam.Gamma = 1.32; // показатель адиабаты c_p/c_v, CH4
    solParam.CFL = 0.9; // число Куранта
    solParam.MaxIter = 1e7; // максимальное кол-во итераций
    solParam.Ma = 0.1; // число Маха

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
    GodunovSolver solver(CH4, solParam, SystemOfEquationType::shockwave1, RiemannSolverType::HLLESolver);
    double h = 1;
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));

    solver.setBorderConditions(&borderConditionShockwave);  // for shockwave

    solver.solve();
    //////////////////////////////////////////////////////////////

}