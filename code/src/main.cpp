/**
 * @file main.cpp
 * @brief Main file for the numerical wave solver.
 *
 * This file contains the main functions and logic for the numerical wave
 * solver. It initializes the model parameters, sets up the grid and time
 * variables for either staggered or collocated grids, and runs the
 * time-stepping loop for the simulation.
 *
 *
 *  that initializes the run parameters using the iput steering file,
 * defines the grid and time variables for either staggered or collocated grids,
 * sets up the flux, dispersion, and wavemaker objects, and runs the
 * time-stepping loop. saves the simulation results to binary files, which will
 * be read and plotted using the python scripts.
 *
 * @details
 * The main components of the code include:
 * - **Input Handling**: Reads user input parameters from the steering file
 * created by the python script.
 * - **Grid Definition**: Sets up the grid based on the selected scheme
 * (staggered or collocated).
 * - **Time Variables**: Allocates and initializes time variables for the
 * simulation. This includes the free surface elevation and velocity.
 * - **Flux Computation**: Computes the SW fluxes based on the selected
 * numerical scheme (HLLC, HLL, central upwind or upwind for the conservative
 * staggered scheme).
 * - **Dispersion**: Implements and discretizes the dispersion terms based on
 * the selected grid type.
 * - **Wavemaker**: Generates wave forcing terms based on the approach described
 * in wei et al. (1999).
 * - **Time Integration**: Solves the governing equations using Runge-Kutta
 * methods (1st, 2nd, or 3rd order).
 * - **Output**: Saves time variables to binary files for post-processing and
 * visualization.
 *
 * @note The code uses OpenMP for parallelization and requires the OpenMP
 * library to be linked during compilation.
 */

#include <algorithm>
#include <cstring>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <vector>

#include "Dispersion.h"
#include "FluxComputation.h"
#include "ReadUserInput.h"
#include "SolveEquation.h"
#include "TimeVariables.h"
#include "Wavemaker.h"

int main(int argc, char **argv) {
  // read steering file
  std::string nametest(argv[1]);
  Input In(nametest);
  In.LOGOUT();

  // define time variables
  State *S0, *S1, *S2, *S3;
  S0 = new State(In);
  S1 = new State(In);
  if (In.OrderTime >= 2)
    S2 = new State(In);
  if (In.OrderTime >= 3)
    S3 = new State(In);

  // define fluxes
  flux *Fcollocated;
  flux_staggered *Fstaggered;

  if (S0->Grid == grid::staggered) {
    Fstaggered = new flux_staggered(S0->N);
  } else {
    Fcollocated = new flux(S0->N);
  }

  // define dispersion
  Dispersion Disp(*S0);

  // define wave maker
  WaveMaker WM(In, S0->N);

  // print info
  printINFO(In);

  // time loop
  float t1 = 0.0f;
  unsigned int K1 = 0;
  float dt1;

  unsigned int Kprint = 0;
  float Tprint = 0;

  unsigned int Kgauges = 0;
  float Tgauges = 0;
  std::vector<float> *Hgauge = new std::vector<float>[In.IndexGauges.size()];

  for (unsigned int i = 0; i < In.IndexGauges.size(); i++) {
    if (In.IndexGauges[i] + positionWM > S0->N) {
      std::cout << "Error : Gauge position out of domain" << '\n';
      exit(1);
    }
    In.IndexGauges[i] = In.IndexGauges[i] + positionWM;
    Hgauge[i].push_back(In.IndexGauges[i]); // save gauges position
  }

  // start computation
  std::cout << std::endl;
  std::cout << "Computation Progress: " << std::endl;

  double start = omp_get_wtime();

  while (t1 <= In.Time) {

    // Compute wave maker source term
    WM.ComputeSourceTerm(t1, In.dx);

    // compute SWE fluxes
    if (S0->Grid == grid::staggered) {
      Fstaggered->ComputeFlux(*S0, In.OrderReconstruction);
    } else {
      Fcollocated->ComputeFlux(*S0, In.Scheme, In.OrderReconstruction);
    }

    // CFL condition for time step
    float maxWS;
    if (S0->Grid == grid::staggered) {
      maxWS = *std::max_element(Fstaggered->WS, Fstaggered->WS + S0->N);
    } else {
      maxWS = *std::max_element(Fcollocated->WS, Fcollocated->WS + S0->N);
    }

    dt1 = In.CourantNumber * S0->dx / maxWS;

    // Runge-Kutta time stepping
    if (S0->Grid == grid::staggered) {
      if (In.OrderTime == 1) {
        SolveEquationStaggered_step1(*S0, *S1, In, *Fstaggered, Disp, WM, dt1);
      }

      else if (In.OrderTime == 2) {
        SolveEquationStaggered_step1(*S0, *S1, In, *Fstaggered, Disp, WM, dt1);
        SolveEquationStaggered_step2(*S0, *S1, *S2, In, *Fstaggered, Disp, WM,
                                     dt1);
      }

      else if (In.OrderTime == 3) {
        SolveEquationStaggered_step1(*S0, *S1, In, *Fstaggered, Disp, WM, dt1);
        SolveEquationStaggered_step2(*S0, *S1, *S2, In, *Fstaggered, Disp, WM,
                                     dt1);
        SolveEquationStaggered_step3(*S0, *S2, *S3, In, *Fstaggered, Disp, WM,
                                     dt1);
      }

      else {
        std::cout << "Error : OrderTime not implemented" << '\n';
        exit(1);
      }
    }

    else {
      if (In.OrderTime == 1) {
        SolveEquationCollocated_step1(*S0, *S1, In, *Fcollocated, Disp, WM,
                                      dt1);
      }

      else if (In.OrderTime == 2) {
        SolveEquationCollocated_step1(*S0, *S1, In, *Fcollocated, Disp, WM,
                                      dt1);
        SolveEquationCollocated_step2(*S0, *S1, *S2, In, *Fcollocated, Disp, WM,
                                      dt1);
      }

      else if (In.OrderTime == 3) {
        SolveEquationCollocated_step1(*S0, *S1, In, *Fcollocated, Disp, WM,
                                      dt1);
        SolveEquationCollocated_step2(*S0, *S1, *S2, In, *Fcollocated, Disp, WM,
                                      dt1);
        SolveEquationCollocated_step3(*S0, *S2, *S3, In, *Fcollocated, Disp, WM,
                                      dt1);
      }

      else {
        std::cout << "Error : OrderTime not implemented" << '\n';
        exit(1);
      }
    }

    // Print variables

    if (In.OrderTime == 1) {
      if (In.IndexGauges.size() != 0)
        RecordTimeSeries(*S0, *S1, In, Hgauge, t1, dt1, Kgauges, Tgauges);
      PrintVariables(*S0, *S1, In, t1, dt1, Kprint, Tprint);
      S0->Update(*S1);
    }

    else if (In.OrderTime == 2) {
      if (In.IndexGauges.size() != 0)
        RecordTimeSeries(*S0, *S2, In, Hgauge, t1, dt1, Kgauges, Tgauges);
      PrintVariables(*S0, *S2, In, t1, dt1, Kprint, Tprint);
      S0->Update(*S2);
    }

    else if (In.OrderTime == 3) {
      if (In.IndexGauges.size() != 0)
        RecordTimeSeries(*S0, *S3, In, Hgauge, t1, dt1, Kgauges, Tgauges);
      PrintVariables(*S0, *S3, In, t1, dt1, Kprint, Tprint);
      S0->Update(*S3);
    }

    t1 += dt1;
    K1++;

    if (t1 >= 0 && t1 <= In.Time * 1.5)
      printProgress(t1 / In.Time);
    else {
      std::cout << "Error : Computational Instability" << '\n';
      exit(1);
    }
  }

  double end = omp_get_wtime();

  // save time series
  if (In.IndexGauges.size() != 0)
    PrintTimeSeries(Hgauge, In);

  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "The numerical simulation took " << end - start << " s"
            << std::endl;

  return 0;
}