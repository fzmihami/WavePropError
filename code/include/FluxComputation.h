#pragma once

#include "TimeVariables.h"

class flux
{
    private:
        unsigned int N;

    public:
        float *WS; // Wave speed for the CFL computation 
        float *Hn, *Qn;
        float *HL, *HR, *QL, *QR;
        float *FluxH, *FluxU;

    public:
        flux(unsigned int N);
        void ComputeFlux(const State& state, const scheme& Scheme, const unsigned int& OrderReconstruction);
        ~flux();
};

class flux_staggered
{
    private:
        unsigned int N;

    public:
        float *WS;  // Wave speed for the CFL computation
        float *HUn; // Hn*Un approximation in the continuity equation
        float *HUUn; // Hn*Un^2 approximation in the momentum equation

    public:
        flux_staggered(unsigned int N);
        void ComputeFlux(const State& state, const unsigned int& OrderReconstruction);
        ~flux_staggered();
    
};

// Flux Reconstruction
void Recontruction_Order1(float* Hn, float* Hl, float *Hr, unsigned int N);
void Recontruction_Order2(float* Hn, float* Hl, float *Hr, unsigned int N);
void Recontruction_Order3(float* Hn, float* Hl, float *Hr, unsigned int N);
void Recontruction_Order5(float* Hn, float* Hl, float *Hr, unsigned int N);
float minmod(float a, float b, float c);
float upwind(float V, float Q1, float Q2);
float upwind(float V, float Q1, float Q2, float Q3, float Q4);

// Flux Computation
void Flux_CentralUpwind(float* Hl, float *Hr, float* Ql, float *Qr, float *FluxH, float *FluxQ, unsigned int N);
void Flux_HLL(float* Hl, float *Hr, float* Ql, float *Qr, float *FluxH, float *FluxQ, unsigned int N);
void Flux_HLLC(float* Hl, float *Hr, float* Ql, float *Qr, float *FluxH, float *FluxQ, unsigned int N);

