#pragma once

#include <iostream>
#include "TimeVariables.h"
#include "GlobalVariables.h"

class Dispersion
{
    private:
        unsigned int N;

    public:
        float *PhiC, *PhiM; // PhiC: Dispersive term for the continuity equation, PhiM: Dispersive term for the momentum equation
        float *Uxx; // Second derivative of the velocity field

        float *D1, *D2, *D3; // Diaginal terms of the matrix P
    
    public:
        Dispersion(const State& S0);
        void ComputeDispersion(State& S0);
        void SolveTridiagonalMatrix(State& S1);
        ~Dispersion();
};

void ThomasAlgorithm(const float* a, const float* b, const float* c, float* d, int n);