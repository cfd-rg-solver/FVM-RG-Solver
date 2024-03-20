//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
#include <filesystem>
#include <energycalc.h>


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

    ///////////////////////////////////////////////////////////////
    ///////////////////// Molecula data //////////////////////////
    ///////////////////////////  CH4  ////////////////////////////
    ///
    // swtting CH4 mixture
    MixtureComponent methane;
    methane.name = "CH4";
    methane.molarMass = 0.016043;
    methane.mass = 2.663732314e-26;
    methane.epsilonDevK = 151.4; // have written in .doc file on drive // epsilon/kB
    methane.sigma = 3.737e-10; // m
    methane.D_diss = 36685.823189; // cm^-1, converted from 438.86 kJ/mol
    methane.numberAtoms = 5;
    methane.numberOfModes = 4;
    methane.omega_eByMode = { 3025.5, 1582.7, 3156.8, 1367.4 }; // sm^-1, all other data, related with length, is in m!
    methane.numberVibrLvlByMode = { 4, 4, 4, 4 }; // { 10, 18, 9, 21 }; 
    methane.dByMode = { 1, 2, 3, 3 };

    for (int i1 = 0; i1 < methane.numberVibrLvlByMode[0]; i1++)
    {
        for (int i2 = 0; i2 < methane.numberVibrLvlByMode[1]; i2++)
        {
            for (int i3 = 0; i3 < methane.numberVibrLvlByMode[2]; i3++)
            {
                for (int i4 = 0; i4 < methane.numberVibrLvlByMode[3]; i4++)
                {
                    double e = (
                        methane.omega_eByMode[0] * (i1 + methane.dByMode[0] / 2.) +
                        methane.omega_eByMode[1] * (i2 + methane.dByMode[1] / 2.) +
                        methane.omega_eByMode[2] * (i3 + methane.dByMode[2] / 2.) +
                        methane.omega_eByMode[3] * (i4 + methane.dByMode[3] / 2.)
                        );
                    if (e < methane.D_diss) {
                        std::vector<int> inds = { i1, i2, i3, i4 };
                        methane.possibleVibrInds.push_back(inds);
                    }
                }
            }
        }
    }

    double viscocity_methane = 1.1123e-05; // for below set sonditions


    std::vector<MixtureComponent> tmp3 = { methane };
    Mixture CH4(tmp3);

    //////////////////////////////////////////////////////////////
    //////////////////// MODEL ///////////////////////////////////
    //////////////////////////////////////////////////////////////
    OneTempApproxMultiModes oneTempApproxMultiModes;

    ///////////////////////////////////////////////////////////////
    ///////////////// Border Condition for Shock Wave ////////////
    ////////////////////////// CH4 ///////////////////////////////
    ///
    // рассматриваем уравнения граничных условий, пусть left = 0, right = n:

    // add calculation based on Mach number

    // Data for simulation:
    // Ma = 3
    // T = 300 K
    // p = 100 Pa
    // (speed of sound 450.06 m/s)
    //
    double velocity_left = 1350; // shock wave, so the velocity is supersonic, let's set it to ??? m/s ~ 5 Ma for methane at room temperature
    double density_left = 4.283556702e-05; // kg/m^3, calculated for atmospheric pressure
    double T_left = 300; // Kelvin
    double pressure_left = UniversalGasConstant * T_left * density_left / methane.molarMass;

    double velocity_right = 308.3333;
    double density_right = 0.0001875;
    double T_right = 688.99176;
    double pressure_right = UniversalGasConstant * T_right * density_right / methane.molarMass;

    BorderConditionShockwave borderConditionShockwave;
    borderConditionShockwave.setBorderParameters(
        velocity_left, density_left, T_left,
        velocity_right, density_right, T_right
    );

    borderConditionShockwave.setEnergyCalculator(&oneTempApproxMultiModes);
    //////////////////////////////////////////////////////////////
    ////////////////// Start param for Shockwave /////////////////
    ////////////////////////////  CH4  ///////////////////////////

    ShockwaveGapDistribution startParamShockwaveCH4;
    macroParam leftStartParam(CH4);
    macroParam rightStartParam(CH4);

    leftStartParam.density = density_left;
    leftStartParam.pressure = pressure_left;
    leftStartParam.velocity_tau = velocity_left;
    leftStartParam.velocity_normal = 0;
    leftStartParam.velocity = velocity_left;
    leftStartParam.fractionArray[0] = 1;
    leftStartParam.densityArray[0] = leftStartParam.fractionArray[0] * leftStartParam.density;

    rightStartParam.density = density_right;
    rightStartParam.pressure = pressure_right;
    rightStartParam.velocity_normal = 0;
    rightStartParam.velocity_tau = velocity_right;
    rightStartParam.velocity = velocity_right;
    rightStartParam.fractionArray[0] = 1;
    rightStartParam.densityArray[0] = rightStartParam.fractionArray[0] * rightStartParam.density;

    startParamShockwaveCH4.setDistributionParameter(leftStartParam, rightStartParam);
    startParamShockwaveCH4.setEnergyCalculator(&oneTempApproxMultiModes);
    //////////////////////////////////////////////////////////////
    ///////////////// Solver param for Shockwave ////////////////
    ////////////////////////////////////////////////////////////

    solverParams solParam;
    solParam.NumCell     = 0.5e2/2 + 2; // Число расчтеных ячеек с учетом двух фиктивных ячеек
    // solParam.Gamma       = 1.67;        // Ar
    // solParam.Gamma       = 1.32;        // O2_O
    solParam.Gamma       = 1.304;       // CH4, but its also implemented changable in macroparam
    solParam.CFL         = 0.9;         // Число Куранта
    solParam.MaxIter     = 5000;        // максимальное кол-во итераций
    solParam.Ma          = 3;         // Число Маха

    double precision = 1E-7; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(10000);


    DataWriter writer(outputData);
    DataReader reader(outputData + "\prev_data");

    reader.read();
    vector<macroParam> startParameters;
    reader.getPoints(startParameters);


    // GodunovSolver solver(Ar ,solParam, SystemOfEquationType::shockwave1, RiemannSolverType::HLLESolver);
    GodunovSolver solver(CH4, solParam, SystemOfEquationType::shockwave2, RiemannSolverType::HLLESolver);

    double M_PI = 3.14159265358979323846;
    // double MFP = viscocity_argon / pressure_left * sqrt(M_PI * UniversalGasConstant * T_left / argon.molarMass); // mean free path length for argon
    double MFP = viscocity_methane / pressure_left * sqrt(M_PI * UniversalGasConstant * T_left / methane.molarMass); // mean free path length for methane
    std::cout << "mean free path: " << MFP << std::endl;
    double h = 30 * MFP; // m
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));

    solver.setEnergyCalculator(&oneTempApproxMultiModes);
    solver.setBorderConditions(&borderConditionShockwave);

    // solver.setStartDistribution(&startParamShockwaveAr);
    solver.setStartDistribution(&startParamShockwaveCH4);

    solver.solve();
}
