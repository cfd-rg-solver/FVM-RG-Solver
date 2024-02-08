#include "godunovsolver.h"
#include "bordercondition.h"
#include "mixture.h"
#include <Python.h>
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
    std::vector<MixtureComponent> tmp1 = { methane }; // - однокомпонентная смесь
    Mixture CH4(tmp1);
    /////////////////////

    //////////////////////////////////////////////////////////////
    ////////////////// Start param for Shockwave /////////////////
    ////////////////////////////  CH4  ///////////////////////////
    UniformDistributionBorder startParamShockwaveCH4;
    UniformDistributionBorderPersonal startParamShockwaveCH4Personal;
    macroParam startParamCH4(CH4);

    startParamCH4.temp = 200; // Kelvin
    startParamCH4.velocity_normal = 0;

    startParamCH4.density = 0.971; // kg/m^3, 200K https://www.engineeringtoolbox.com/methane-density-specific-weight-temperature-pressure-d_2020.html
    startParamCH4.fractionArray[0] = 1; // однокомпонентная смесь
    startParamCH4.densityArray[0] = startParamCH4.fractionArray[0] * startParamCH4.density;

    startParamShockwaveCH4.setBorderCondition(&borderConditionShockwave);
    startParamShockwaveCH4.setDistributionParameter(startParamCH4);
    //////////////////////////////////////////////////////////////

    
    //////////////////////////////////////////////////////////////
    ///////////////// Border Condition for Shockwave ////////////
    //////////////////////////// Ar ////////////////////////////

    /// Ar (argon) /// 
    MixtureComponent argon;
    argon.name = "argon";
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.numberAtoms = 1;
    // argon.numberVibrLvl = ?;
    // argon.epsilonDevK = ?
    // argon.sigma = 3.33E-10; ?
    // argon.omega_e = ?
    /////////////////////

    // рассматриваем уравнения граничных условий,
    // пусть left = 0, right = n: 
    double velocity_left = 0;
    double density_left = 0.800773; // todo какое ???
    double T_left = 1000;
    double pressure_left = UniversalGasConstant * T_left * density_left / argon.molarMass; // todo норм ???
    double energy_left = 3 * kB * T_left / (2 * argon.mass); // одноатомный аргон

    Py_Initialize();
    PyObject *name, *load_module, *func, *callfunc, *args;
    name = PyUnicode_FromString((char*)"nonlinear_border_system_solver");
    load_module = PyImport_Import(name);
    func = PyObject_GetAttrString(load_module, (char*)"solver");
    args = PyTuple_Pack(4,
        PyFloat_FromDouble(velocity_left),
        PyFloat_FromDouble(density_left),
        PyFloat_FromDouble(pressure_left),
        PyFloat_FromDouble(energy_left)
    );
    callfunc = PyObject(func, args);
    vector<double> solver_out = PyTuple_AsVector(vallfunc);
    Py_Finalize();

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
    UniformDistributionBorder startParamShockwaveAr;
    UniformDistributionBorderPersonal startParamShockwaveArPersonal;
    macroParam startParamAr(Ar);

    startParamAr.temp = 200; // Kelvin
    startParamAr.velocity_normal = 0;

    startParamAr.density = 2.409; // kg/m^3, 200K https://www.engineeringtoolbox.com/argon-density-specific-weight-temperature-pressure-d_2089.html
    startParamAr.fractionArray[0] = 1; // однокомпонентная смесь
    startParamCH4.densityArray[0] = startParamAr.fractionArray[0] * startParamAr.density;

    startParamShockwaveAr.setBorderCondition(&borderConditionShockwave);
    startParamShockwaveAr.setDistributionParameter(startParamAr);
    //////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////
    /////////////////////// General Solver ///////////////////////
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
    // GodunovSolver solver(CH4, solParam, SystemOfEquationType::shockwave2, RiemannSolverType::HLLESolver);
    GodunovSolver solver(Ar, solParam, SystemOfEquationType::shockwave1, RiemannSolverType::HLLESolver);
    double h = 1;
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));
    solver.setBorderConditions(&borderConditionShockwave);  // for shockwave
    solver.solve();
    //////////////////////////////////////////////////////////////

}