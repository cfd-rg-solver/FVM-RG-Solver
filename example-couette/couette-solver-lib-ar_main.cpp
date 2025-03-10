//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
#include <filesystem>


std::string GetCurrentWorkingDir( void ) {
    std::filesystem::path currentWorkingDir = std::filesystem::current_path();
    std::filesystem::path parentDir = currentWorkingDir.parent_path().parent_path();

   std::string res = parentDir.string() + "/FVM-RG-Solver/example-couette";

    return res;
}

namespace fs = std::filesystem;
int main()
{

    std::string outputData = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << outputData << std::endl;

    //////////////////////////////////////////////////////////////
    ///////////////////// Border Condition for Couette ///////////
    //////////////////////////////////////////////////////////////
    int caseType = 6;
    double T_up_wall;
    double T_down_wall;
    double velocity_up;
    double velocity_down;

    double argonSpeedOfSound = 307.73; // m/s	
    double aronSpeedOfSound2 = 295.10; // m/s at 251 K

    if(caseType == 0)
    {
        T_up_wall = 273;
        T_down_wall = 273;
        velocity_up = 0; //300;
        velocity_down = 0;
    }
    else if(caseType == 1)
    {
        T_up_wall = 273;
        T_down_wall = 273;
        velocity_up = 1888.84;
        velocity_down = 0;
    }

    else if(caseType == 2)
    {
        T_up_wall = 273;
        T_down_wall = 273;
        velocity_up = 1 * argonSpeedOfSound;
        velocity_down = 0;
    }

    else if(caseType == 3)
    {
        T_up_wall = 273;
        T_down_wall = 273;
        velocity_up = 2 * argonSpeedOfSound;
        velocity_down = 0;
    }

    else if(caseType == 4)
    {
        T_up_wall = 273;
        T_down_wall = 273;
        velocity_up = 3 * argonSpeedOfSound;
        velocity_down = 0;
    }

    else if(caseType == 5)
    {
        T_up_wall = 273;
        T_down_wall = 273;
        velocity_up = 4.5 * argonSpeedOfSound;
        velocity_down = 0;
    }

    else if(caseType == 6)
    {
        T_up_wall = 251.05;
        T_down_wall = 251.05;
        velocity_up = 1556.16;
        velocity_down = 0;
    }


    BorderConditionCouette borderConditionCouette;
    borderConditionCouette.setWallParameters(velocity_up, velocity_down, T_up_wall, T_down_wall);

    //////////////////////////////////////////////////////////////
    ///////////////////// Border Condition for Couette ///////////
    /////////////////////////// Slip /////////////////////////////

    BorderConditionCouetteSlip borderConditionCouetteSlip;
    borderConditionCouetteSlip.setWallParameters(velocity_up, velocity_down, T_up_wall, T_down_wall);
    double accommodationCoeff = 1.0;
    borderConditionCouetteSlip.setAccommodationCoeff(accommodationCoeff);

    //////////////////////////////////////////////////////////////

    // Ar
    MixtureComponent argon;
    argon.name = "Ar";
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21/kB;
    argon.numberAtoms = 1;
    argon.sigma = 3.33E-10;


    std::vector<MixtureComponent> tmp = {argon};
    Mixture Ar(tmp);

    //////////////////////////////////////////////////////////////
    ///////////////////// Start param for Couette ////////////////
    ////////////////////////////  Ar  ///////////////////////////

    UniformDistributionBorder startParamCouetteAr;
    UniformDistributionBorder startParamCouetteArSlip; // Slip Border
    startParamCouetteAr.setMixture(Ar); // TODO temp
    startParamCouetteArSlip.setMixture(Ar); // TODO temp
    macroParam startParamAr(Ar);
    int newSolving = 2;

    double pressure;

    if(newSolving == 1)
    {
        startParamAr.density = 0.000115; // 0.0000115; //  correct case
        startParamAr.fractionArray[0] = 1;
        startParamAr.densityArray[0] =  startParamAr.fractionArray[0] * startParamAr.density;

        startParamAr.temp = 273; 
        startParamAr.velocity_tau = 0;
        startParamAr.velocity_normal = 0;

        pressure = startParamAr.density * T_up_wall * UniversalGasConstant / argon.molarMass;

        startParamCouetteAr.setBorderCondition(&borderConditionCouette);
        startParamCouetteAr.setDistributionParameter(startParamAr);

        startParamCouetteArSlip.setBorderCondition(&borderConditionCouetteSlip);
        startParamCouetteArSlip.setDistributionParameter(startParamAr);

        std::cout << "pressure: " << pressure << std::endl;
    }
    else if (newSolving == 0)
    {
        DataWriter writer(outputData); 
        DataReader reader(outputData + "/prev-data");

        reader.read();
        vector<macroParam> startParameters;
        reader.getPoints(startParameters);

        startParamCouetteAr.setBorderCondition(&borderConditionCouette);
        startParamCouetteAr.setDistributionParameter(startParameters);

        startParamCouetteArSlip.setBorderCondition(&borderConditionCouetteSlip);
        startParamCouetteArSlip.setDistributionParameter(startParameters);

        pressure = startParameters[0].density * T_up_wall * UniversalGasConstant / argon.molarMass; // Pa, for the above set of conditions

        std::cout << "pressure: " << pressure << std::endl;
    }
    else if (newSolving == 2)
    {
        startParamAr.density = 0.003852;
        startParamAr.fractionArray[0] = 1;
        startParamAr.densityArray[0] =  startParamAr.fractionArray[0] * startParamAr.density;

        startParamAr.temp = 251; 
        startParamAr.velocity_tau = 0;
        startParamAr.velocity_normal = 0;

        pressure = startParamAr.density * T_up_wall * UniversalGasConstant / argon.molarMass;

        startParamCouetteAr.setBorderCondition(&borderConditionCouette);
        startParamCouetteAr.setDistributionParameter(startParamAr);

        startParamCouetteArSlip.setBorderCondition(&borderConditionCouetteSlip);
        startParamCouetteArSlip.setDistributionParameter(startParamAr);

        std::cout << "pressure: " << pressure << std::endl;
    }

    //////////////////////////////////////////////////////////////

    solverParams solParam;
    solParam.NumCell     = 202;    // Число расчтеных ячеек с учетом двух фиктивных ячеек
    solParam.Gamma    = 1.67;   // Ar
    solParam.CFL      = 0.95;    // Число Куранта 0.9
    solParam.MaxIter     = 8000000; // максимальное кол-во итареций
    solParam.Ma       = 0.1;    // Число маха

    double precision = 1E-7; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(10000);


    // DataWriter writer(outputData);
    DataWriter writer(outputData);
    DataReader reader(outputData + "/prev_data");

    reader.read();
    vector<macroParam> startParameters;
    reader.getPoints(startParameters);


    double viscocity_argon = 2.0988e-05;
    double viscocity_argon2 = 1.9532e-05;

    GodunovSolver solver(Ar, solParam, SystemOfEquationType::couette1, RiemannSolverType::HLLESolver);
    
    double MFP = viscocity_argon2 / pressure * sqrt(M_PI * UniversalGasConstant * T_up_wall / 2. / argon.molarMass); // Mean free path length for Argon

    double h = 0.1; // General length between the walls
    // h = h - 4 * MFP; // Effective length between the walls (length without Knudsen layers)

    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));

    

    bool BCSlip = 0;
    if(BCSlip)
    {
        solver.setBorderConditions(&borderConditionCouetteSlip); // Slip border
        solver.setStartDistribution(&startParamCouetteArSlip); // Slip border
        // writer.writeSimulationParam("border type", "slip");
        // writer.writeSimulationParam("accommodationCoeff",accommodationCoeff);
    }
    else
    {
        solver.setBorderConditions(&borderConditionCouette);
        solver.setStartDistribution(&startParamCouetteAr);
        // writer.writeSimulationParam("border type", "noSlip");
    }

    // writer.writeSimulationParam("density", startParamAr.density);
    // writer.writeSimulationParam("number of cell", solParam.NumCell);
    // writer.writeSimulationParam("gamma", solParam.Gamma);
    // writer.writeSimulationParam("CFL", solParam.CFL);

    std::cout << "mean free path: " << MFP << std::endl;

    std::cout << "Knudsen number: " << MFP / h << std::endl;

    std::cout << "cell size: " << h/(solParam.NumCell - 2) << std::endl;

    solver.solve();
}
