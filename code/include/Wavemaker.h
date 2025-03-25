# pragma once

#include <iostream>

#include "GlobalVariables.h"
#include "ReadUserInput.h"

class WaveMaker
{
    public:
        unsigned int Nf; // Number of frequencies
        unsigned int N; // Number of grid points
        float Beta; 
        unsigned int posWM; // Position of the wavemaker
        float depth; // Water depth
        float width; // half the width of the wavemaker
        unsigned int Ns; // Number of grid points for the wavemaker

        float Kp; // peak wavelength

    private:
        const float de = 0.5f;
        const float alpha3 = alpha1 + float(1.0/3.0);

    public:
        float *Freq, *Amp, *Phase;
        float *Omg, *K;
        float *Di;

        float *SourceTerm;

    public:
        WaveMaker(const Input& In, const unsigned int& N);
        void ComputeSourceTerm(const float& t, const float& dx);
        ~WaveMaker();
};

float SecondorderEquation(float a, float b, float c);
float wavelength(float omega, float depth);