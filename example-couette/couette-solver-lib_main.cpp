//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
#include <filesystem>


std::string GetCurrentWorkingDir( void ) {
    std::filesystem::path currentWorkingDir = std::filesystem::current_path();
    std::filesystem::path parentDir = currentWorkingDir.parent_path().parent_path();

    std::string res = parentDir.string() + "/FVM-RG-Solver/main/example-couette";

//    std::string res = parentDir.string() + "/main/example-couette"; // !normal case using qt
    return res;
//    return currentWorkingDir; // without qt
//    return currentWorkingDir.string();
}

//#define Liia
namespace fs = std::filesystem;
int main()
{

<<<<<<< HEAD:example/couette-solver-lib_main.cpp
    // Ar
=======
    std::string outputData = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << outputData << std::endl;

>>>>>>> origin/main:example-couette/couette-solver-lib_main.cpp
    MixtureComponent argon;
    argon.name = "Ar";
    argon.density = 1.7839e-3; // case 1
    //argon.density = 0.800773;
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21/kB; //! Mistake
    argon.sigma = 3.33E-10;
    Mixture Ar({argon});

    // O2/O

    MixtureComponent Oxygen2;
    Oxygen2.name = "O2";
    Oxygen2.density = 1.42897e-3;
    Oxygen2.molarMass = 0.032;
    Oxygen2.mass = 5.313725014956E-26;
    Oxygen2.epsilonDevK = 107.4;
    Oxygen2.sigma = 3.458E-10;

    MixtureComponent Oxygen;
    Oxygen.name = "O";
    Oxygen.density = 1.42897e-3;
    Oxygen.molarMass = 0.015999;
    Oxygen.mass = 5.313725014956E-26 / 2.;
    Oxygen.epsilonDevK = 80;
    Oxygen.sigma = 2.75E-10;

    std::vector<MixtureComponent> tmp = {Oxygen2,Oxygen};
    Mixture O2_O(tmp);

    /////Выбор смеси/////

    Mixture mixture = O2_O;

    /////////////////////

    macroParam startParam(mixture);
    startParam.fractionArray[0] = 0.5;
    startParam.densityArray[0] =  startParam.fractionArray[0] * mixture.components[0].density;
    if(mixture.NumberOfComponents == 2)
    {
        startParam.fractionArray[1] = 0.5;
        startParam.densityArray[1] =  startParam.fractionArray[1] * mixture.components[1].density;
    }
    startParam.density = 0;
    for(int i = 0; i < mixture.NumberOfComponents; i++)
        startParam.density += startParam.densityArray[i];
    startParam.temp = 1000; //140
    startParam.velocity_tau = 0.00001;


    solverParams solParam;
    solParam.NumCell     = 102;    // Число расчтеных ячеек с учетом двух фиктивных ячеек
    solParam.Gamma    = 1.67;    // Показатель адиабаты
    solParam.CFL      = 0.9;    // Число Куранта
    solParam.MaxIter     = 10000000; // максимальное кол-во итареций
    solParam.Ma       = 0.1;    // Число маха

    double precision = 1E-5; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(10000);

<<<<<<< HEAD:example/couette-solver-lib_main.cpp
    // это меняешь под себя. Он так создаст папку data
    // если не использовать setWriter, то записи не будет, но папка создастся, ибо она в конструкторе зашита
    // он автоматически очищает папку перед новым рассчётом
#ifdef Liia
    DataWriter writer("D:/main/work_materials/CdExamples/couette-mark/couette-solver-lib");
    DataReader reader("D:/main/work_materials/CdExamples/couette-mark/couette-solver-lib/prev_data");
#else
    DataWriter writer("D:/study/course_work");
    DataReader reader("D:/study/course_work/prev_data");
#endif
=======
    DataWriter writer(outputData);
    DataReader reader(outputData + "\prev_data");

>>>>>>> origin/main:example-couette/couette-solver-lib_main.cpp
    reader.read();
    vector<macroParam> startParameters;
    reader.getPoints(startParameters);
    double T1wall = 1000;
    double T2wall = 1000;
    double velocity = 300;
    double h = 1;
    GodunovSolver solver(mixture,startParam,solParam, SystemOfEquationType::couette2Alt, RiemannSolverType::HLLESolver);
    writer.setDelta_h(h / (solParam.NumCell - 2));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell - 2));
    solver.setBorderConditions(velocity,T2wall,T1wall);
    solver.solve();
}
