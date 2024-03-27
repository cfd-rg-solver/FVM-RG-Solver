//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
#include <filesystem>
#include <energycalc.h>


std::string GetCurrentWorkingDir(void) {
    std::filesystem::path currentWorkingDir = std::filesystem::current_path();
    std::filesystem::path parentDir = currentWorkingDir.parent_path().parent_path();

    std::string res = parentDir.string() + "/FVM-RG-Solver/example-shockwave"; // ! win case

    // std::string res = parentDir.string() + "/main/example-couette"; // !normal case using qt
    return res;
    //    return currentWorkingDir; // without qt
    //    return currentWorkingDir.string();
}

double M_PI = 3.14159265358979323846;

namespace fs = std::filesystem;
int main()
{

    std::string outputData = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << outputData << std::endl;

    ///////////////////////////////////////////////////////////////
    ///////////////////// Molecula data //////////////////////////
    ///////////////////////////  Ar  ////////////////////////////
    ///
    // swtting Ar mixture
    MixtureComponent argon;
    argon.name = "Ar";
    argon.molarMass = 0.039948; // kg/mol
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21 / kB; //! Mistake
    argon.sigma = 3.33E-10;
    argon.numberAtoms = 1;
    argon.omega_e = 0;

    std::vector<MixtureComponent> tmp2 = { argon };
    Mixture Ar(tmp2);

    double viscocity_argon = 2.2725e-05; // for below set sonditions

    ///////////////////////////////////////////////////////////////
    ///////////////////// Molecula data //////////////////////////
    ///////////////////////////  CH4  ////////////////////////////
    ///
    // swtting CH4 mixture
    /*
    MixtureComponent methane;
    methane.name = "CH4";
    methane.molarMass = 0.016043;
    methane.mass = 2.663732314e-26;
    methane.epsilonDevK = 151.4; // have written in .doc file on drive // epsilon/kB
    methane.sigma = 3.737e-10; // m
    methane.D_diss = 3668582.3189; // m^-1!, converted from 438.86 kJ/mol
    methane.numberAtoms = 5;
    methane.numberOfModes = 4;
    methane.omega_eByMode = { 302550, 158270, 315680, 136740 }; // m^-1! all other data, related with length, is in m!
    methane.numberVibrLvlByMode = { 1,1,1,1 };//{ 10, 18, 9, 21 }; //{ 4, 4, 4, 4 };
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
    */

    //////////////////////////////////////////////////////////////
    //////////////////// MODEL ///////////////////////////////////
    //////////////////////////////////////////////////////////////
    // OneTempApproxMultiModes oneTempApproxMultiModes; // CH4
    OneTempApprox oneTempApprox; // Ar

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
    // METHANE SET:
    /*
    double velocity_left = 1350.18;
    double density_left = 0.00064318; // kg/m^3, calculated for atmospheric pressure
    double T_left = 300; // Kelvin
    double pressure_left = UniversalGasConstant * T_left * density_left / methane.molarMass;

    double velocity_right = 359.36775;
    double density_right = 0.002416;
    double T_right = 766.889; // ! from solver (methane case)
    double pressure_right = UniversalGasConstant * T_right * density_right / methane.molarMass;
    */
    // ARGON SET:
    double velocity_left = 1225.88;
    double density_left = 0.000106663; // kg/m^3 ??? calculated as pressure*molarMass / (R * T_left), but its has different dimention in NIST
    double T_left = 300; // Kelvin
    double pressure_left = UniversalGasConstant * T_left * density_left / argon.molarMass;

    double velocity_right = 370.2480;
    double density_right = 0.000353158;
    double T_right = 1612.9346; // ! from approx solver, other one is not fixed for methane case
    double pressure_right = UniversalGasConstant * T_right * density_right / argon.molarMass;

    BorderConditionShockwave borderConditionShockwave;
    borderConditionShockwave.setBorderParameters(
        velocity_left, density_left, T_left,
        velocity_right, density_right, T_right
    );

    // borderConditionShockwave.setEnergyCalculator(&oneTempApproxMultiModes); // CH4
    borderConditionShockwave.setEnergyCalculator(&oneTempApprox); // Ar

    //////////////////////////////////////////////////////////////
    ////////////////// Start param for Shockwave /////////////////
    ////////////////////////////  CH4  ///////////////////////////
    /*
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
    */
    //////////////////////////////////////////////////////////////
    ////////////////// Start param for Shockwave /////////////////
    ////////////////////////////  Ar  ///////////////////////////
    ShockwaveGapDistribution startParamShockwaveAr;
    macroParam leftStartParam(Ar);
    macroParam rightStartParam(Ar);

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

    startParamShockwaveAr.setDistributionParameter(leftStartParam, rightStartParam);
    startParamShockwaveAr.setEnergyCalculator(&oneTempApprox);

    //////////////////////////////////////////////////////////////
    ///////////////// Solver param for Shockwave ////////////////
    ////////////////////////////////////////////////////////////

    solverParams solParam;
    solParam.NumCell = 30 + 2; // Число расчтеных ячеек с учетом двух фиктивных ячеек
    solParam.Gamma = 1.67;        // Ar
    // solParam.Gamma       = 1.32;        // O2_O
    // solParam.Gamma       = 1.304;       // CH4, but its also implemented changable in macroparam
    solParam.CFL = 0.9;         // Число Куранта
    solParam.MaxIter = 10000;        // максимальное кол-во итераций
    solParam.Ma = 3;         // Число Маха, сейчас не влияет на решатель, просто формальность

    double precision = 1E-7; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(10000);


    DataWriter writer(outputData);
    DataReader reader(outputData + "\prev_data");

    reader.read();
    vector<macroParam> startParameters;
    reader.getPoints(startParameters);


    GodunovSolver solver(Ar, solParam, SystemOfEquationType::shockwave1, RiemannSolverType::HLLESolver);
    //GodunovSolver solver(CH4, solParam, SystemOfEquationType::shockwave2, RiemannSolverType::HLLESolver);

    double MFP = viscocity_argon / pressure_left * sqrt(M_PI * UniversalGasConstant * T_left / argon.molarMass); // mean free path length for argon
    //double MFP = viscocity_methane / pressure_left * sqrt(M_PI * UniversalGasConstant * T_left / methane.molarMass); // mean free path length for methane
    std::cout << "mean free path: " << MFP << std::endl;
    double h = 40 * MFP; // m
    std::cout << "considering h = MFP * " << h / MFP << std::endl;
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));

    solver.setEnergyCalculator(&oneTempApprox); // Ar
    //solver.setEnergyCalculator(&oneTempApproxMultiModes); // CH4
    solver.setBorderConditions(&borderConditionShockwave);

    solver.setStartDistribution(&startParamShockwaveAr);
    // solver.setStartDistribution(&startParamShockwaveCH4);

    solver.solve();
}
