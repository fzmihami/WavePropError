/**
 * @file SolveEquation.cpp
 * @brief Implementation of functions for solving the conservative form of the
 * Nwogu's equation.
 *
 * This file contains the implementation of functions for solving the Nwogu's
 * equations using both staggered and collocated grids. The functions are
 * designed to handle different time-stepping methods and numerical schemes.
 */

#include <cmath>
#include <iostream>
#include <string.h>

#include "Boundaries.h"
#include "SolveEquation.h"

void SolveEquationStaggered_step1(State &S0, State &S1, const Input &In,
                                  flux_staggered &Fstaggered, Dispersion &Disp,
                                  WaveMaker &WM, const float &dt1) {
  // Compute dispersion
  Disp.ComputeDispersion(S0);

// Solve Continuity equation
#pragma omp parallel for
  for (unsigned int i = 0; i < S0.N; i++) {
    S1.Hn[i] = S0.Hn[i] -
               dt1 * (Fstaggered.HUn[i + 1] - Fstaggered.HUn[i]) * S0.dx_1 -
               dt1 * Disp.PhiC[i] + dt1 * WM.SourceTerm[i];
  }

// Compute new Eta
#pragma omp parallel for
  for (unsigned int i = 0; i < S0.N; i++) {
    S1.Eta[i] = S1.Hn[i] + S1.Zn[i];
  }

  // Solve Momentum equation
  float coeff = In.OrderTime == 1 ? 1.0 : 0.0;
#pragma omp parallel for
  for (unsigned int i = 1; i < S0.N; i++) {
    float H1 = 0.5 * (S1.Hn[i] + S1.Hn[i - 1]);
    float H0 = 0.5 * (S0.Hn[i] + S0.Hn[i - 1]);

    S1.Pn[i] =
        S0.Pn[i] -
        dt1 * (Fstaggered.HUUn[i] - Fstaggered.HUUn[i - 1]) * S0.dx_1 -
        dt1 * S0.Un[i] * 0.5 * (Disp.PhiC[i] + Disp.PhiC[i - 1]) +
        dt1 * (H1 - H0) * dt1 * Disp.PhiM[i] -
        coeff * dt1 * G * S1.Hn[i] * (S1.Eta[i] - S1.Eta[i - 1]) * S0.dx_1;
  }
  S1.Pn[0] = S1.Pn[1];
  S1.Pn[S0.N] = S1.Pn[S0.N - 1];

  // Compute new Un
  Disp.SolveTridiagonalMatrix(S1);

  // Apply sponge layer
  ApplySpongeLayer(S1);
}

void SolveEquationStaggered_step2(State &S0, State &S1, State &S2,
                                  const Input &In, flux_staggered &Fstaggered,
                                  Dispersion &Disp, WaveMaker &WM,
                                  const float &dt1) {
  float RKcoeff1 = In.OrderTime == 2 ? 0.5 : 3.0f / 4.0f;
  float RKcoeff2 = In.OrderTime == 2 ? 0.5 : 1.0f / 4.0f;

  // Compute fluxes
  Fstaggered.ComputeFlux(S1, In.OrderReconstruction);

  // Compute dispersion
  Disp.ComputeDispersion(S1);

// Solve Continuity equation
#pragma omp parallel for
  for (unsigned int i = 0; i < S0.N; i++) {
    S2.Hn[i] = RKcoeff1 * S0.Hn[i] + RKcoeff2 * S1.Hn[i] +
               RKcoeff2 * (-dt1 * (Fstaggered.HUn[i + 1] - Fstaggered.HUn[i]) *
                               S0.dx_1 -
                           dt1 * Disp.PhiC[i] + dt1 * WM.SourceTerm[i]);
  }

// Compute new Eta
#pragma omp parallel for
  for (unsigned int i = 0; i < S0.N; i++) {
    S2.Eta[i] = S2.Hn[i] + S2.Zn[i];
  }

  // Solve Momentum equation
  float coeff = In.OrderTime == 2 ? 1.0 : 0.0;
#pragma omp parallel for
  for (unsigned int i = 1; i < S0.N; i++) {
    float H1 = 0.5 * (S2.Hn[i] + S2.Hn[i - 1]);
    float H0 = 0.5 * (S1.Hn[i] + S1.Hn[i - 1]);

    S2.Pn[i] =
        RKcoeff1 * S0.Pn[i] + RKcoeff2 * S1.Pn[i] +
        RKcoeff2 *
            (-dt1 * (Fstaggered.HUUn[i] - Fstaggered.HUUn[i - 1]) * S0.dx_1 -
             dt1 * S1.Un[i] * 0.5 * (Disp.PhiC[i] + Disp.PhiC[i - 1]) +
             dt1 * (H1 - H0) * dt1 * Disp.PhiM[i]) -
        coeff * dt1 * G * S2.Hn[i] * (S2.Eta[i] - S2.Eta[i - 1]) * S0.dx_1;
  }
  S2.Pn[0] = S2.Pn[1];
  S2.Pn[S0.N] = S2.Pn[S0.N - 1];

  // Compute new Un
  Disp.SolveTridiagonalMatrix(S2);

  // Apply sponge layer
  ApplySpongeLayer(S2);
}

void SolveEquationStaggered_step3(State &S0, State &S2, State &S3,
                                  const Input &In, flux_staggered &Fstaggered,
                                  Dispersion &Disp, WaveMaker &WM,
                                  const float &dt1) {
  float RKcoeff1 = 1.0f / 3.0f;
  float RKcoeff2 = 2.0f / 3.0f;

  // Compute fluxes
  Fstaggered.ComputeFlux(S2, In.OrderReconstruction);

  // Compute dispersion
  Disp.ComputeDispersion(S2);

// Solve Continuity equation
#pragma omp parallel for
  for (unsigned int i = 0; i < S0.N; i++) {
    S3.Hn[i] = RKcoeff1 * S0.Hn[i] + RKcoeff2 * S2.Hn[i] +
               RKcoeff2 * (-dt1 * (Fstaggered.HUn[i + 1] - Fstaggered.HUn[i]) *
                               S0.dx_1 -
                           dt1 * Disp.PhiC[i] + dt1 * WM.SourceTerm[i]);
  }

// Compute new Eta
#pragma omp parallel for
  for (unsigned int i = 0; i < S0.N; i++) {
    S3.Eta[i] = S3.Hn[i] + S3.Zn[i];
  }

// Solve Momentum equation
#pragma omp parallel for
  for (unsigned int i = 1; i < S0.N; i++) {
    float H1 = 0.5 * (S3.Hn[i] + S3.Hn[i - 1]);
    float H0 = 0.5 * (S2.Hn[i] + S2.Hn[i - 1]);

    S3.Pn[i] =
        RKcoeff1 * S0.Pn[i] + RKcoeff2 * S2.Pn[i] +
        RKcoeff2 *
            (-dt1 * (Fstaggered.HUUn[i] - Fstaggered.HUUn[i - 1]) * S0.dx_1 -
             dt1 * S2.Un[i] * 0.5 * (Disp.PhiC[i] + Disp.PhiC[i - 1]) +
             dt1 * (H1 - H0) * dt1 * Disp.PhiM[i]) -
        dt1 * G * S3.Hn[i] * (S3.Eta[i] - S3.Eta[i - 1]) * S0.dx_1;
  }
  S3.Pn[0] = S3.Pn[1];
  S3.Pn[S0.N] = S3.Pn[S0.N - 1];

  // Compute new Un
  Disp.SolveTridiagonalMatrix(S3);

  // Apply sponge layer
  ApplySpongeLayer(S3);
}

// Collocated grid /////////////////////////////////////////////////////////////

void SolveEquationCollocated_step1(State &S0, State &S1, const Input &In,
                                   flux &Fcollocated, Dispersion &Disp,
                                   WaveMaker &WM, const float &dt1) {
  // Compute dispersion
  Disp.ComputeDispersion(S0);

// Solve Continuity equation
#pragma omp parallel for
  for (unsigned int i = 1; i < S0.N - 1; i++) {
    S1.Hn[i] =
        S0.Hn[i] -
        dt1 * (Fcollocated.FluxH[i] - Fcollocated.FluxH[i - 1]) * S0.dx_1 -
        dt1 * Disp.PhiC[i] + dt1 * WM.SourceTerm[i];
  }
  S1.Hn[0] = S1.Hn[1];
  S1.Hn[S0.N - 1] = S1.Hn[S0.N - 2];

// Compute new Eta
#pragma omp parallel for
  for (unsigned int i = 0; i < S0.N; i++) {
    S1.Eta[i] = S1.Hn[i] + S1.Zn[i];
  }

// Solve Momentum equation
#pragma omp parallel for
  for (unsigned int i = 1; i < S0.N - 1; i++) {
    S1.Pn[i] =
        S0.Pn[i] -
        dt1 * (Fcollocated.FluxU[i] - Fcollocated.FluxU[i - 1]) * S0.dx_1 -
        dt1 * S0.Un[i] * Disp.PhiC[i] +
        dt1 * (S1.Hn[i] - S0.Hn[i]) * dt1 * Disp.PhiM[i];
  }
  S1.Pn[0] = S1.Pn[1];
  S1.Pn[S0.N] = S1.Pn[S0.N - 1];

  // Compute new Un
  Disp.SolveTridiagonalMatrix(S1);

  // Apply sponge layer
  ApplySpongeLayer(S1);
}

void SolveEquationCollocated_step2(State &S0, State &S1, State &S2,
                                   const Input &In, flux &Fcollocated,
                                   Dispersion &Disp, WaveMaker &WM,
                                   const float &dt1) {
  float RKcoeff1 = In.OrderTime == 2 ? 0.5 : 3.0f / 4.0f;
  float RKcoeff2 = In.OrderTime == 2 ? 0.5 : 1.0f / 4.0f;

  // Compute fluxes
  Fcollocated.ComputeFlux(S1, In.Scheme, In.OrderReconstruction);

  // Compute dispersion
  Disp.ComputeDispersion(S1);

// Solve Continuity equation
#pragma omp parallel for
  for (unsigned int i = 1; i < S0.N - 1; i++) {
    S2.Hn[i] =
        RKcoeff1 * S0.Hn[i] + RKcoeff2 * S1.Hn[i] +
        RKcoeff2 * (-dt1 * (Fcollocated.FluxH[i] - Fcollocated.FluxH[i - 1]) *
                        S0.dx_1 -
                    dt1 * Disp.PhiC[i] + dt1 * WM.SourceTerm[i]);
  }
  S2.Hn[0] = S2.Hn[1];
  S2.Hn[S0.N - 1] = S2.Hn[S0.N - 2];

// Compute new Eta
#pragma omp parallel for
  for (unsigned int i = 0; i < S0.N; i++) {
    S2.Eta[i] = S2.Hn[i] + S2.Zn[i];
  }

// Solve Momentum equation
#pragma omp parallel for
  for (unsigned int i = 1; i < S0.N - 1; i++) {
    S2.Pn[i] =
        RKcoeff1 * S0.Pn[i] + RKcoeff2 * S1.Pn[i] +
        RKcoeff2 * (-dt1 * (Fcollocated.FluxU[i] - Fcollocated.FluxU[i - 1]) *
                        S0.dx_1 -
                    dt1 * S1.Un[i] * Disp.PhiC[i] +
                    dt1 * (S2.Hn[i] - S1.Hn[i]) * dt1 * Disp.PhiM[i]);
  }
  S2.Pn[0] = S2.Pn[1];
  S2.Pn[S0.N] = S2.Pn[S0.N - 1];

  // Compute new Un
  Disp.SolveTridiagonalMatrix(S2);

  // Apply sponge layer
  ApplySpongeLayer(S2);
}

void SolveEquationCollocated_step3(State &S0, State &S2, State &S3,
                                   const Input &In, flux &Fcollocated,
                                   Dispersion &Disp, WaveMaker &WM,
                                   const float &dt1) {
  float RKcoeff1 = 1.0f / 3.0f;
  float RKcoeff2 = 2.0f / 3.0f;

  // Compute fluxes
  Fcollocated.ComputeFlux(S2, In.Scheme, In.OrderReconstruction);

  // Compute dispersion
  Disp.ComputeDispersion(S2);

// Solve Continuity equation
#pragma omp parallel for
  for (unsigned int i = 1; i < S0.N - 1; i++) {
    S3.Hn[i] =
        RKcoeff1 * S0.Hn[i] + RKcoeff2 * S2.Hn[i] +
        RKcoeff2 * (-dt1 * (Fcollocated.FluxH[i] - Fcollocated.FluxH[i - 1]) *
                        S0.dx_1 -
                    dt1 * Disp.PhiC[i] + dt1 * WM.SourceTerm[i]);
  }
  S3.Hn[0] = S3.Hn[1];
  S3.Hn[S0.N - 1] = S3.Hn[S0.N - 2];

// Compute new Eta
#pragma omp parallel for
  for (unsigned int i = 0; i < S0.N; i++) {
    S3.Eta[i] = S3.Hn[i] + S3.Zn[i];
  }

// Solve Momentum equation
#pragma omp parallel for
  for (unsigned int i = 1; i < S0.N - 1; i++) {
    S3.Pn[i] =
        RKcoeff1 * S0.Pn[i] + RKcoeff2 * S2.Pn[i] +
        RKcoeff2 * (-dt1 * (Fcollocated.FluxU[i] - Fcollocated.FluxU[i - 1]) *
                        S0.dx_1 -
                    dt1 * S2.Un[i] * Disp.PhiC[i] +
                    dt1 * (S3.Hn[i] - S2.Hn[i]) * dt1 * Disp.PhiM[i]);
  }
  S3.Pn[0] = S3.Pn[1];
  S3.Pn[S0.N] = S3.Pn[S0.N - 1];

  // Compute new Un
  Disp.SolveTridiagonalMatrix(S3);

  // Apply sponge layer
  ApplySpongeLayer(S3);
}