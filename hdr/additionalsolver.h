#pragma once

#include <vector>
#include <mutex>

#include "global.h"
#include <utility>

typedef std::pair<vector<Matrix>, vector<Matrix>> conservative;
struct AdditionalSolver
{
public:
    macroParam ExacRiemanSolver(macroParam left, macroParam right, double Gamma);
    conservative SolveEvolutionExplFirstOrder(vector<Matrix> F1, Matrix F2, Matrix F3,  vector<Matrix> U1old, Matrix U2old, Matrix U3old, double dt, double delta_h);

};
