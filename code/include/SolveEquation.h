#pragma once

#include "ReadUserInput.h"
#include "TimeVariables.h"
#include "FluxComputation.h"
#include "Dispersion.h"
#include "Wavemaker.h"


void SolveEquationStaggered_step1(State& S0, State& S1, const Input& In, flux_staggered& Fstaggered, Dispersion& Disp, WaveMaker& WM, const float& dt1);

void SolveEquationStaggered_step2(State& S0, State& S1, State& S2, const Input& In, flux_staggered& Fstaggered, Dispersion& Disp, WaveMaker& WM, const float& dt1);

void SolveEquationStaggered_step3(State& S0, State& S2, State& S3, const Input& In, flux_staggered& Fstaggered, Dispersion& Disp, WaveMaker& WM, const float& dt1);


void SolveEquationCollocated_step1(State& S0, State& S1, const Input& In, flux& Fcollocated, Dispersion& Disp, WaveMaker& WM, const float& dt1);

void SolveEquationCollocated_step2(State& S0, State& S1, State& S2, const Input& In, flux& Fcollocated, Dispersion& Disp, WaveMaker& WM, const float& dt1);

void SolveEquationCollocated_step3(State& S0, State& S2, State& S3, const Input& In, flux& Fcollocated, Dispersion& Disp, WaveMaker& WM, const float& dt1);
